"""CLI command for building cluster DB from dat-bypass pipeline output.

Maps cluster assignments back to original spectra via scan_titles.txt files
generated during parquet_to_dat conversion. Uses DuckDB for efficient joins
and aggregations.
"""
import gc
import re
import click
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pyspectrafuse.common.constant import ParquetSchemaAdapter
from pyspectrafuse.common.schemas import CLUSTER_META_SCHEMA, PSM_MEMBERSHIP_SCHEMA
from pyspectrafuse.common.duckdb_ops import (
    compute_stats_and_purity,
    scan_parquet_scalars,
    scan_parquet_arrays,
    write_parquet_from_df,
    _detect_scan_expr,
)
from pyspectrafuse.commands.spectrum2msp import find_target_ext_files

logger = logging.getLogger(__name__)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def parse_cluster_tsv(tsv_path: str) -> pd.DataFrame:
    """Parse MaRaCluster cluster TSV."""
    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=['mgf_path', 'scannr', 'cluster_id'],
                     skip_blank_lines=False)
    df = df.dropna()
    df['scannr'] = df['scannr'].astype(int)
    df['cluster_id'] = df['cluster_id'].astype(str)
    # Strip window suffix (e.g., _w0, _w1) from MGF basename so it matches
    # the original scan_titles filename (without window index)
    df['mgf_basename'] = df['mgf_path'].apply(
        lambda p: re.sub(r'_w\d+\.mgf$', '.mgf', Path(p).name))
    return df[['mgf_basename', 'scannr', 'cluster_id']]


def parse_scan_titles(titles_path: str, dataset_name: str) -> pd.DataFrame:
    """Parse scan_titles.txt → DataFrame mapping scannr to original spectrum identity."""
    rows = []
    title_pattern = re.compile(r'id=mzspec::(.+?):scan:(\d+):(.+?)$')

    with open(titles_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t', 2)
            if len(parts) < 3:
                continue
            scannr = int(parts[1])
            title = parts[2]

            m = title_pattern.match(title)
            if m:
                run_file = m.group(1)
                orig_scan = int(m.group(2))
                pep_charge = m.group(3)
                last_slash = pep_charge.rfind('/')
                if last_slash > 0:
                    peptidoform = pep_charge[:last_slash]
                    charge = int(pep_charge[last_slash + 1:])
                else:
                    peptidoform = pep_charge
                    charge = 0
                rows.append((scannr, run_file, orig_scan, peptidoform, charge, dataset_name))

    return pd.DataFrame(rows, columns=['scannr', 'run_file_name', 'orig_scan',
                                        'peptidoform', 'charge', 'dataset'])


def build_cluster_db(
    cluster_tsv: str,
    scan_titles_files: List[str],
    parquet_dirs: List[str],
    dataset_names: List[str],
    species: str,
    instrument: str,
    charge_str: str,
    output_dir: str,
    method_type: str = 'best',
    chunk_size: int = 200_000,
) -> Tuple[Optional[str], Optional[str]]:
    """Build cluster DB from dat-bypass output using scan_titles mapping.

    Uses DuckDB for efficient parquet scanning and aggregation:
      Pass 1: DuckDB SEMI JOIN to read scalars (PEP, q-value, precursor_mz)
      Pass 2: DuckDB SEMI JOIN to read spectrum arrays for best PSM per cluster
    """
    charge_int = int(charge_str.replace('charge', ''))
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Parse cluster TSV ──
    cluster_df = parse_cluster_tsv(cluster_tsv)
    logger.info(f"Cluster TSV: {len(cluster_df)} entries, "
                f"{cluster_df['cluster_id'].nunique()} clusters")

    # ── Parse scan_titles and build membership ──
    titles_dfs = []
    for titles_file, ds_name in zip(scan_titles_files, dataset_names):
        titles_df = parse_scan_titles(titles_file, ds_name)
        stem = Path(titles_file).name.replace('.scan_titles.txt', '')
        titles_df['mgf_basename'] = f'{stem}.mgf'
        titles_dfs.append(titles_df)

    all_titles = pd.concat(titles_dfs, ignore_index=True)
    membership = cluster_df.merge(all_titles, on=['mgf_basename', 'scannr'], how='inner')

    if len(membership) < len(cluster_df):
        logger.warning(f"{len(cluster_df) - len(membership)} spectra lost in join")
    logger.info(f"Membership: {len(membership)} PSMs in {membership['cluster_id'].nunique()} clusters")

    # ── Write membership to temp parquet for DuckDB-only pipeline ──
    # Instead of loading everything into pandas, write membership as parquet
    # and let DuckDB handle the enrichment + PSM writing entirely in SQL.
    import tempfile
    tmp_dir = tempfile.mkdtemp(prefix='cluster_db_')
    membership_tmp = str(Path(tmp_dir) / 'membership.parquet')
    membership.to_parquet(membership_tmp, index=False)
    del membership; gc.collect()

    # Collect all parquet file paths
    all_parquet_files = []
    for pdir in parquet_dirs:
        all_parquet_files.extend(find_target_ext_files(pdir, '.parquet'))

    # ── Pass 1: Build PSM membership entirely in DuckDB ──
    # DuckDB reads source parquet with predicate pushdown (charge filter)
    # and SEMI JOINs with the membership temp parquet — no pandas intermediary.
    logger.info(f"Pass 1: Building PSM membership via DuckDB (charge={charge_int})...")

    # Use a file-backed DuckDB database so large JOINs can spill to disk
    # instead of OOM-killing the process.
    duckdb_path = str(Path(tmp_dir) / 'build.duckdb')
    conn = duckdb.connect(duckdb_path)
    try:
        conn.execute("SET memory_limit='4GB'")
        conn.execute(f"SET temp_directory='{tmp_dir}'")
        conn.execute("SET threads=4")
    except Exception:
        conn.close()
        raise

    # Build the source parquet UNION with column normalization
    source_parts = []
    for ppath in all_parquet_files:
        cols_df = conn.execute(f"SELECT name FROM parquet_schema('{ppath}')").fetchdf()
        available = set(cols_df['name'])

        run_col = 'run_file_name' if 'run_file_name' in available else 'reference_file_name'
        charge_col = 'charge' if 'charge' in available else 'precursor_charge'

        mz_parts = []
        if 'observed_mz' in available: mz_parts.append('observed_mz')
        if 'calculated_mz' in available: mz_parts.append('calculated_mz')
        if 'exp_mass_to_charge' in available: mz_parts.append('exp_mass_to_charge')
        mz_expr = f"COALESCE({', '.join(mz_parts)}, 0.0)" if mz_parts else '0.0'

        pep_col = 'posterior_error_probability' if 'posterior_error_probability' in available else 'NULL'

        # Handle scan column type
        scan_expr = _detect_scan_expr(conn, ppath)

        if 'global_qvalue' in available:
            qvalue_expr = 'global_qvalue'
        elif 'additional_scores' in available:
            qvalue_expr = """(
                SELECT s.score_value
                FROM UNNEST(additional_scores) AS t(s)
                WHERE s.score_name = 'global_qvalue'
                LIMIT 1
            )"""
        else:
            qvalue_expr = 'NULL'

        source_parts.append(f"""
            SELECT
                {run_col} AS run_file_name,
                {scan_expr} AS orig_scan,
                {mz_expr} AS precursor_mz,
                {pep_col} AS posterior_error_probability,
                {qvalue_expr} AS global_qvalue
            FROM read_parquet('{ppath}')
            WHERE {charge_col} = {charge_int}
        """)

    source_union = " UNION ALL ".join(source_parts)

    # Write PSM membership directly to parquet via DuckDB
    psm_path = str(out_dir / 'psm_cluster_membership.parquet')
    try:
        conn.execute(f"""
            COPY (
                SELECT
                    m.cluster_id,
                    'mzspec:' || m.dataset || ':' || m.run_file_name || ':scan:'
                        || CAST(m.orig_scan AS VARCHAR) || ':' || m.peptidoform
                        || '/' || CAST({charge_int} AS VARCHAR) AS usi,
                    m.dataset AS project_accession,
                    m.run_file_name AS reference_file_name,
                    CAST(m.orig_scan AS INTEGER) AS scan,
                    m.peptidoform,
                    CAST({charge_int} AS TINYINT) AS charge,
                    COALESCE(s.precursor_mz, 0.0) AS precursor_mz,
                    s.posterior_error_probability,
                    s.global_qvalue,
                    '{species}' AS species,
                    '{instrument}' AS instrument
                FROM read_parquet('{membership_tmp}') m
                LEFT JOIN ({source_union}) s
                    ON m.run_file_name = s.run_file_name
                    AND m.orig_scan = s.orig_scan
                ORDER BY m.cluster_id
            ) TO '{psm_path}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """)

        n_psms = conn.execute(f"SELECT COUNT(*) FROM read_parquet('{psm_path}')").fetchone()[0]
        logger.info(f"PSM membership: {n_psms} rows written to {psm_path}")

        # ── Cluster Stats + Purity via DuckDB ──
        stats_purity = compute_stats_and_purity(psm_path)

        # ── Best PSM per cluster + spectrum arrays (chunked to avoid OOM) ──
        logger.info(f"Pass 2: Finding best PSM per cluster...")

        # Get best PSM keys via DuckDB on the membership parquet we just wrote
        best_keys = conn.execute(f"""
            SELECT cluster_id, reference_file_name AS run_file_name,
                   scan AS orig_scan, peptidoform, precursor_mz
            FROM (
                SELECT *,
                       ROW_NUMBER() OVER (
                           PARTITION BY cluster_id
                           ORDER BY posterior_error_probability ASC NULLS LAST
                       ) AS _rn
                FROM read_parquet('{psm_path}')
            )
            WHERE _rn = 1
        """).fetchdf()
        n_clusters = len(best_keys)
        best_keys['orig_scan'] = best_keys['orig_scan'].astype(int)
        logger.info(f"Best PSMs: {n_clusters} clusters")
    finally:
        conn.close()

    # Merge with stats/purity early (small — just scalars)
    best_keys = best_keys.merge(stats_purity, on='cluster_id', how='left')
    best_keys['purity'] = best_keys['purity'].fillna(1.0).astype(np.float32)
    best_keys['member_count'] = best_keys['member_count'].fillna(1).astype(np.int32)
    best_keys['project_count'] = best_keys['project_count'].fillna(1).astype(np.int16)
    del stats_purity; gc.collect()

    # ── Write cluster metadata in chunks (avoids loading all arrays at once) ──
    meta_path = str(out_dir / 'cluster_metadata.parquet')
    writer = None
    total_written = 0

    def to_f32_list(arr):
        if arr is None:
            return []
        if isinstance(arr, np.ndarray):
            return arr.astype(np.float32).tolist()
        if isinstance(arr, (list, tuple)):
            return [float(x) for x in arr]
        return []

    for chunk_start in range(0, n_clusters, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_clusters)
        chunk_keys = best_keys.iloc[chunk_start:chunk_end]
        logger.info(f"Pass 2 chunk {chunk_start//chunk_size + 1}: "
                     f"clusters {chunk_start}-{chunk_end} of {n_clusters}")

        # Read arrays for this chunk only
        keys_for_join = chunk_keys[['run_file_name', 'orig_scan']].copy()
        arrays_df = scan_parquet_arrays(
            all_parquet_files, keys_for_join, charge_int,
            temp_directory=tmp_dir)

        chunk_with_arrays = chunk_keys.merge(
            arrays_df, on=['run_file_name', 'orig_scan'], how='left')
        del arrays_df, keys_for_join; gc.collect()

        has_mz = 'mz_array' in chunk_with_arrays.columns
        has_int = 'intensity_array' in chunk_with_arrays.columns

        chunk_meta = pd.DataFrame({
            'cluster_id': chunk_with_arrays['cluster_id'].values,
            'species': species,
            'instrument': instrument,
            'charge': np.int8(charge_int),
            'peptidoform': chunk_with_arrays['peptidoform'].values,
            'peptide_sequence': chunk_with_arrays['peptidoform'].str.split('/').str[0].str.replace(
                r'\[.*?\]', '', regex=True).values,
            'consensus_mz_array': chunk_with_arrays['mz_array'].apply(to_f32_list).values if has_mz else [[] for _ in range(len(chunk_with_arrays))],
            'consensus_intensity_array': chunk_with_arrays['intensity_array'].apply(to_f32_list).values if has_int else [[] for _ in range(len(chunk_with_arrays))],
            'consensus_method': method_type,
            'precursor_mz': pd.to_numeric(chunk_with_arrays['precursor_mz'], errors='coerce').fillna(0).values,
            'member_count': chunk_with_arrays['member_count'].values,
            'project_count': chunk_with_arrays['project_count'].values,
            'best_pep': chunk_with_arrays['best_pep'].values if 'best_pep' in chunk_with_arrays.columns else np.nan,
            'best_qvalue': chunk_with_arrays['best_qvalue'].values if 'best_qvalue' in chunk_with_arrays.columns else np.nan,
            'purity': chunk_with_arrays['purity'].values,
        })
        del chunk_with_arrays; gc.collect()

        chunk_table = pa.Table.from_pandas(chunk_meta, schema=CLUSTER_META_SCHEMA, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(meta_path, schema=CLUSTER_META_SCHEMA, compression='zstd')
        writer.write_table(chunk_table)
        total_written += len(chunk_meta)
        del chunk_meta, chunk_table; gc.collect()

    if writer is not None:
        writer.close()

    # Clean up temp files
    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)

    logger.info(f"Cluster metadata: {total_written} clusters, "
                f"{Path(meta_path).stat().st_size / 1e6:.1f} MB")

    return meta_path, psm_path


@click.command("build-cluster-db", short_help="Build cluster DB from dat-bypass output")
@click.option('--cluster_tsv', required=True, type=click.Path(exists=True),
              help='MaRaCluster cluster TSV file')
@click.option('--scan_titles', required=True, multiple=True,
              type=click.Path(exists=True),
              help='scan_titles.txt file(s) from parquet_to_dat (one per dataset)')
@click.option('--parquet_dir', required=True, multiple=True,
              type=click.Path(exists=True),
              help='Parquet directory/directories containing PSM data')
@click.option('--dataset_name', required=True, multiple=True,
              help='Dataset name(s) matching scan_titles order')
@click.option('--species', required=True, help='Species name')
@click.option('--instrument', required=True, help='Instrument name')
@click.option('--charge', required=True, help='Charge value (e.g., charge2)')
@click.option('--output_dir', required=True, type=click.Path(),
              help='Output directory for cluster DB files')
@click.option('--method_type', default='best',
              type=click.Choice(['best', 'bin']),
              help='Consensus spectrum method (default: best)')
@click.option('--chunk_size', default=200_000, type=int,
              help='Clusters per chunk for memory-efficient array loading')
def build_cluster_db_cmd(cluster_tsv, scan_titles, parquet_dir, dataset_name,
                         species, instrument, charge, output_dir,
                         method_type, chunk_size):
    """Build cluster_metadata.parquet and psm_cluster_membership.parquet
    from dat-bypass pipeline output.

    Uses scan_titles.txt files to map MaRaCluster dat indices back to
    original spectra in the parquet files. DuckDB-powered for efficient
    joins and aggregations at any scale.
    """
    if len(scan_titles) != len(dataset_name):
        raise click.ClickException(
            f"Number of --scan_titles ({len(scan_titles)}) must match "
            f"--dataset_name ({len(dataset_name)})")

    meta_path, psm_path = build_cluster_db(
        cluster_tsv=cluster_tsv,
        scan_titles_files=list(scan_titles),
        parquet_dirs=list(parquet_dir),
        dataset_names=list(dataset_name),
        species=species,
        instrument=instrument,
        charge_str=charge,
        output_dir=output_dir,
        method_type=method_type,
        chunk_size=chunk_size,
    )

    if meta_path:
        click.echo(f"Cluster metadata: {meta_path}")
        click.echo(f"PSM membership: {psm_path}")
    else:
        click.echo("No output generated (no matching spectra)")
