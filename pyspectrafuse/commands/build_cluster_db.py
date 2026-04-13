"""CLI command for building cluster DB from dat-bypass pipeline output.

Maps cluster assignments back to original spectra via scan_titles.txt files
generated during parquet_to_dat conversion. Handles the naming mismatch between
dat bypass dummy MGFs and the standard inject_cluster_info() convention.
"""
import gc
import re
import click
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pyspectrafuse.common.constant import ParquetSchemaAdapter
from pyspectrafuse.commands.spectrum2msp import find_target_ext_files

logger = logging.getLogger(__name__)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

# Output schemas
CLUSTER_META_SCHEMA = pa.schema([
    ('cluster_id', pa.string()),
    ('species', pa.string()),
    ('instrument', pa.string()),
    ('charge', pa.int8()),
    ('peptidoform', pa.string()),
    ('peptide_sequence', pa.string()),
    ('consensus_mz_array', pa.list_(pa.float32())),
    ('consensus_intensity_array', pa.list_(pa.float32())),
    ('consensus_method', pa.string()),
    ('precursor_mz', pa.float64()),
    ('member_count', pa.int32()),
    ('project_count', pa.int16()),
    ('best_pep', pa.float64()),
    ('best_qvalue', pa.float64()),
    ('purity', pa.float32()),
])

PSM_MEMBERSHIP_SCHEMA = pa.schema([
    ('cluster_id', pa.string()),
    ('usi', pa.string()),
    ('project_accession', pa.string()),
    ('reference_file_name', pa.string()),
    ('scan', pa.int32()),
    ('peptidoform', pa.string()),
    ('charge', pa.int8()),
    ('precursor_mz', pa.float64()),
    ('posterior_error_probability', pa.float64()),
    ('global_qvalue', pa.float64()),
    ('species', pa.string()),
    ('instrument', pa.string()),
])


def parse_cluster_tsv(tsv_path: str) -> pd.DataFrame:
    """Parse MaRaCluster cluster TSV."""
    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=['mgf_path', 'scannr', 'cluster_id'],
                     skip_blank_lines=False)
    df = df.dropna()
    df['scannr'] = df['scannr'].astype(int)
    df['cluster_id'] = df['cluster_id'].astype(int).astype(str)
    df['mgf_basename'] = df['mgf_path'].apply(lambda p: Path(p).name)
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


def _scan_parquet_scalars(
    needed: Set[Tuple[str, int]],
    parquet_dirs: List[str],
    charge_int: int,
) -> pd.DataFrame:
    """Read scalar columns (PEP, q-value, precursor_mz) from parquet files."""
    read_cols = ['charge', 'observed_mz', 'calculated_mz', 'run_file_name',
                 'scan', 'posterior_error_probability', 'additional_scores']
    parts = []

    for pdir in parquet_dirs:
        parquet_files = find_target_ext_files(pdir, '.parquet')
        for pf_path in parquet_files:
            pf = pq.ParquetFile(pf_path)
            available = set(pf.schema_arrow.names)
            cols = [c for c in read_cols if c in available]

            for batch in pf.iter_batches(batch_size=200_000, columns=cols):
                df = batch.to_pandas()
                df = ParquetSchemaAdapter.adapt(df)

                if 'charge' in df.columns:
                    df = df[df['charge'] == charge_int]
                    if df.empty:
                        continue

                scan_col = df['scan']
                if scan_col.dtype == object:
                    df['scan_int'] = scan_col.apply(
                        lambda s: int(s[0]) if isinstance(s, (list, np.ndarray)) and len(s) > 0
                        else (int(s) if s is not None else 0))
                else:
                    df['scan_int'] = scan_col.fillna(0).astype(int)

                keys = list(zip(df['reference_file_name'].values, df['scan_int'].values))
                mask = pd.array([k in needed for k in keys], dtype=bool)
                df = df[mask]
                if df.empty:
                    continue

                if 'pepmass' in df.columns:
                    prec_mz = pd.to_numeric(df['pepmass'], errors='coerce')
                    if 'calculated_mz' in df.columns:
                        prec_mz = prec_mz.fillna(pd.to_numeric(df['calculated_mz'], errors='coerce'))
                elif 'calculated_mz' in df.columns:
                    prec_mz = pd.to_numeric(df['calculated_mz'], errors='coerce')
                else:
                    prec_mz = pd.Series(0.0, index=df.index)

                parts.append(pd.DataFrame({
                    'run_file_name': df['reference_file_name'].values,
                    'orig_scan': df['scan_int'].values,
                    'precursor_mz': prec_mz.fillna(0.0).values,
                    'posterior_error_probability': pd.to_numeric(
                        df.get('posterior_error_probability', pd.Series(dtype=float)),
                        errors='coerce').reindex(df.index).values,
                    'global_qvalue': pd.to_numeric(
                        df.get('global_qvalue', pd.Series(dtype=float)),
                        errors='coerce').reindex(df.index).values,
                }))

    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


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

    Two-pass memory-efficient approach:
      Pass 1: Read scalars (PEP, q-value, precursor_mz) for all spectra
      Pass 2: Read spectrum arrays in chunks for best PSM per cluster
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
        # mgf_basename matches the dat filename with .mgf extension
        stem = Path(titles_file).name.replace('.scan_titles.txt', '')
        titles_df['mgf_basename'] = f'{stem}.mgf'
        titles_dfs.append(titles_df)

    all_titles = pd.concat(titles_dfs, ignore_index=True)
    membership = cluster_df.merge(all_titles, on=['mgf_basename', 'scannr'], how='inner')

    if len(membership) < len(cluster_df):
        logger.warning(f"{len(cluster_df) - len(membership)} spectra lost in join")
    logger.info(f"Membership: {len(membership)} PSMs in {membership['cluster_id'].nunique()} clusters")

    # ── Pass 1: Scalars ──
    needed = set(zip(membership['run_file_name'].values,
                     membership['orig_scan'].astype(int).values))
    logger.info(f"Pass 1: Reading scalars for {len(needed)} spectra...")

    scalars_df = _scan_parquet_scalars(needed, parquet_dirs, charge_int)
    if scalars_df.empty:
        logger.warning("No parquet data matched!")
        return None, None
    logger.info(f"Pass 1: Found {len(scalars_df)} spectra")

    enriched = membership.merge(scalars_df, on=['run_file_name', 'orig_scan'],
                                how='left', suffixes=('', '_pq'))
    del scalars_df; gc.collect()

    # ── Write PSM Membership ──
    psm_df = pd.DataFrame({
        'cluster_id': enriched['cluster_id'].values,
        'usi': ('mzspec:' + enriched['dataset'].astype(str) + ':' +
                enriched['run_file_name'].astype(str) + ':scan:' +
                enriched['orig_scan'].astype(str) + ':' +
                enriched['peptidoform'].astype(str) + '/' +
                enriched['charge'].astype(str)),
        'project_accession': enriched['dataset'].values,
        'reference_file_name': enriched['run_file_name'].values,
        'scan': enriched['orig_scan'].values.astype(np.int32),
        'peptidoform': enriched['peptidoform'].values,
        'charge': np.full(len(enriched), charge_int, dtype=np.int8),
        'precursor_mz': pd.to_numeric(enriched['precursor_mz'], errors='coerce').fillna(0).values,
        'posterior_error_probability': pd.to_numeric(
            enriched['posterior_error_probability'], errors='coerce').values,
        'global_qvalue': pd.to_numeric(
            enriched['global_qvalue'], errors='coerce').values,
        'species': species,
        'instrument': instrument,
    })

    psm_path = str(out_dir / 'psm_cluster_membership.parquet')
    psm_df = psm_df.sort_values('cluster_id')
    psm_table = pa.Table.from_pandas(psm_df, schema=PSM_MEMBERSHIP_SCHEMA, preserve_index=False)
    pq.write_table(psm_table, psm_path, compression='zstd')
    logger.info(f"PSM membership: {len(psm_df)} rows, {Path(psm_path).stat().st_size / 1e6:.1f} MB")
    del psm_table; gc.collect()

    # ── Cluster Stats ──
    cluster_stats = psm_df.groupby('cluster_id').agg(
        member_count=('usi', 'count'),
        project_count=('project_accession', 'nunique'),
        best_pep=('posterior_error_probability', 'min'),
        best_qvalue=('global_qvalue', 'min'),
    ).reset_index()

    pair_counts = psm_df.groupby(['cluster_id', 'peptidoform']).size().reset_index(name='_cnt')
    mode_per_cluster = pair_counts.groupby('cluster_id')['_cnt'].max()
    total_per_cluster = psm_df.groupby('cluster_id').size()
    purity_series = mode_per_cluster / total_per_cluster
    del psm_df, pair_counts, mode_per_cluster, total_per_cluster; gc.collect()

    # ── Best PSM per cluster ──
    enriched['_pep'] = pd.to_numeric(enriched['posterior_error_probability'], errors='coerce')
    best_idx = enriched.groupby('cluster_id')['_pep'].idxmin()
    best_rows = enriched.loc[best_idx, ['cluster_id', 'run_file_name', 'orig_scan',
                                         'peptidoform', 'precursor_mz']].copy()
    n_clusters = len(best_rows)
    logger.info(f"Best PSMs: {n_clusters} clusters")
    del enriched; gc.collect()

    # ── Pass 2: Chunked array reading + incremental metadata write ──
    meta_path = str(out_dir / 'cluster_metadata.parquet')
    writer = None

    def to_f32_list(arr):
        if arr is None:
            return []
        if isinstance(arr, np.ndarray):
            return arr.astype(np.float32).tolist()
        if isinstance(arr, (list, tuple)):
            return [float(x) for x in arr]
        return []

    n_chunks = max((n_clusters + chunk_size - 1) // chunk_size, 1)
    best_chunks = np.array_split(best_rows, n_chunks)

    read_cols = ['mz_array', 'intensity_array', 'charge', 'run_file_name', 'scan']

    for chunk_i, chunk_df in enumerate(best_chunks):
        chunk_keys = set(zip(chunk_df['run_file_name'].values,
                             chunk_df['orig_scan'].astype(int).values))
        logger.info(f"Pass 2 chunk {chunk_i+1}/{n_chunks}: {len(chunk_keys)} spectra")

        arrays_dict = {}
        for pdir in parquet_dirs:
            parquet_files = find_target_ext_files(pdir, '.parquet')
            for pf_path in parquet_files:
                pf = pq.ParquetFile(pf_path)
                available = set(pf.schema_arrow.names)
                cols = [c for c in read_cols if c in available]

                for batch in pf.iter_batches(batch_size=200_000, columns=cols):
                    df = batch.to_pandas()
                    df = ParquetSchemaAdapter.adapt(df)
                    if 'charge' in df.columns:
                        df = df[df['charge'] == charge_int]
                        if df.empty:
                            continue

                    scan_col = df['scan']
                    if scan_col.dtype == object:
                        df['scan_int'] = scan_col.apply(
                            lambda s: int(s[0]) if isinstance(s, (list, np.ndarray)) and len(s) > 0
                            else (int(s) if s is not None else 0))
                    else:
                        df['scan_int'] = scan_col.fillna(0).astype(int)

                    keys = list(zip(df['reference_file_name'].values, df['scan_int'].values))
                    mask = [k in chunk_keys for k in keys]
                    df = df[mask]
                    if df.empty:
                        continue

                    for _, row in df.iterrows():
                        key = (row['reference_file_name'], row['scan_int'])
                        if key not in arrays_dict:
                            arrays_dict[key] = (row['mz_array'], row['intensity_array'])

        mz_arrays = []
        int_arrays = []
        for _, row in chunk_df.iterrows():
            key = (row['run_file_name'], int(row['orig_scan']))
            if key in arrays_dict:
                mz_arr, int_arr = arrays_dict[key]
                mz_arrays.append(to_f32_list(mz_arr))
                int_arrays.append(to_f32_list(int_arr))
            else:
                mz_arrays.append([])
                int_arrays.append([])
        del arrays_dict; gc.collect()

        chunk_meta = pd.DataFrame({
            'cluster_id': chunk_df['cluster_id'].values,
            'species': species,
            'instrument': instrument,
            'charge': np.int8(charge_int),
            'peptidoform': chunk_df['peptidoform'].values,
            'peptide_sequence': chunk_df['peptidoform'].str.split('/').str[0].str.replace(
                r'\[.*?\]', '', regex=True).values,
            'consensus_mz_array': mz_arrays,
            'consensus_intensity_array': int_arrays,
            'consensus_method': method_type,
            'precursor_mz': pd.to_numeric(chunk_df['precursor_mz'], errors='coerce').fillna(0).values,
        })
        del mz_arrays, int_arrays

        chunk_meta = chunk_meta.merge(cluster_stats, on='cluster_id', how='left')
        chunk_meta['purity'] = chunk_meta['cluster_id'].map(purity_series).fillna(1.0).astype(np.float32)
        chunk_meta['member_count'] = chunk_meta['member_count'].fillna(1).astype(np.int32)
        chunk_meta['project_count'] = chunk_meta['project_count'].fillna(1).astype(np.int16)

        chunk_table = pa.Table.from_pandas(chunk_meta, schema=CLUSTER_META_SCHEMA, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(meta_path, CLUSTER_META_SCHEMA, compression='zstd')
        writer.write_table(chunk_table)
        del chunk_meta, chunk_table; gc.collect()

    if writer is not None:
        writer.close()

    del best_rows, cluster_stats, purity_series; gc.collect()
    logger.info(f"Cluster metadata: {n_clusters} clusters, "
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
    original spectra in the parquet files. Memory-efficient chunked approach
    for partitions with millions of clusters.
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
