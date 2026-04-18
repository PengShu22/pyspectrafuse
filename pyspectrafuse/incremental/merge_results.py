"""Merge incremental clustering results into existing parquet tables.

When ``build-cluster-db`` is invoked with ``--existing_metadata`` /
``--existing_membership``, this module handles the merge: appending new PSMs
to membership (dedup by USI), updating consensus spectra for existing clusters
that received new members, building metadata for brand-new clusters, and
accumulating provenance.
"""
import glob
import logging
import re
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pyspectrafuse.common.schemas import CLUSTER_META_SCHEMA
from pyspectrafuse.common.duckdb_ops import (
    append_membership_dedup,
    compute_stats_and_purity,
)

logger = logging.getLogger(__name__)


def append_psm_membership(
    existing_path: str,
    new_psms_df: pd.DataFrame,
    output_path: Optional[str] = None,
) -> str:
    """Append new PSM rows to the membership parquet, deduplicating by USI."""
    if output_path is None:
        output_path = existing_path
    if 'resolved_cluster_id' in new_psms_df.columns:
        new_psms_df = new_psms_df.rename(columns={'resolved_cluster_id': 'cluster_id'})
    return append_membership_dedup(existing_path, new_psms_df, output_path)


def _aggregate_source_datasets(membership_path: str) -> pd.DataFrame:
    """Aggregate DISTINCT project_accession per cluster from the membership parquet."""
    import duckdb
    conn = duckdb.connect(':memory:')
    try:
        return conn.execute(f"""
            SELECT cluster_id,
                   list_distinct(array_agg(project_accession)) AS source_datasets
            FROM read_parquet('{membership_path}')
            GROUP BY cluster_id
        """).fetchdf()
    finally:
        conn.close()


def merge_into_existing_db(
    cluster_tsv: str,
    scan_titles_files: List[str],
    existing_metadata_path: str,
    existing_membership_path: str,
    new_parquet_paths: List[str],
    species: str,
    instrument: str,
    charge: str,
    method_type: str,
    output_dir: str,
    project_accession: str,
) -> tuple:
    """Merge a new round of clustering into an existing cluster DB.

    Workflow:
    1. Parse scan_titles; split reps (``rep:{cluster_id}``) from new-data titles.
    2. Resolve new MaRaCluster cluster IDs → either an old cluster_id (if a rep
       lands in the new cluster) or a fresh UUID.
    3. Join new-data titles with the new project's parquet via DuckDB.
    4. Append new PSMs to psm_cluster_membership.parquet (dedup by USI).
    5. For existing clusters receiving new members: swap in better-qvalue
       spectrum (BEST) or recompute consensus (BIN).
    6. For brand-new clusters: build metadata rows from scratch.
    7. Stamp provenance columns: ``is_reused_cluster`` and ``source_datasets``.

    Returns:
        Tuple of (out_metadata_path, out_membership_path).
    """
    from pyspectrafuse.incremental.resolve_clusters_dat import (
        resolve_incremental_clusters,
        load_new_spectra_with_arrays,
    )
    from pyspectrafuse.commands.spectrum2msp import create_consensus_strategy

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    out_membership = str(Path(output_dir) / 'psm_cluster_membership.parquet')
    out_metadata = str(Path(output_dir) / 'cluster_metadata.parquet')

    # Step 1+2: resolve cluster IDs
    id_map, new_spectra_df = resolve_incremental_clusters(cluster_tsv, scan_titles_files)
    logger.info(f"Resolved {len(id_map)} cluster IDs, {len(new_spectra_df)} new spectra")

    # Step 3: enrich new spectra with parquet PSM data
    charge_int = int(str(charge).replace('charge', ''))
    new_df = load_new_spectra_with_arrays(
        new_rows=new_spectra_df,
        parquet_paths=new_parquet_paths,
        charge_int=charge_int,
        project_accession=project_accession,
    )
    if new_df.empty:
        logger.warning("No new spectra matched cluster assignments")
        return out_metadata, out_membership

    # Step 4: append to PSM membership
    psm_new = pd.DataFrame({
        'cluster_id': new_df['cluster_accession'].astype(str),
        'usi': new_df['usi'].astype(str),
        'project_accession': project_accession,
        'reference_file_name': new_df['reference_file_name'].astype(str),
        'scan': new_df['scan'].astype(int),
        'peptidoform': new_df['peptidoform'].astype(str),
        'charge': new_df['charge'].astype(int),
        'precursor_mz': pd.to_numeric(new_df['pepmass'], errors='coerce').astype(float),
        'posterior_error_probability': pd.to_numeric(
            new_df['posterior_error_probability'], errors='coerce').astype(float),
        'global_qvalue': pd.to_numeric(new_df['global_qvalue'], errors='coerce').astype(float),
        'species': species,
        'instrument': instrument,
    })
    append_psm_membership(existing_membership_path, psm_new, out_membership)

    # Step 5: recompute stats + purity and update metadata
    stats_purity = compute_stats_and_purity(out_membership)
    existing_meta = pd.read_parquet(existing_metadata_path).copy()
    existing_meta['cluster_id'] = existing_meta['cluster_id'].astype(str)
    existing_ids = set(existing_meta['cluster_id'])

    new_df_copy = new_df.copy()
    new_df_copy['cluster_accession'] = new_df_copy['cluster_accession'].astype(str)

    new_for_existing_df = new_df_copy[new_df_copy['cluster_accession'].isin(existing_ids)]
    new_brand_new_df = new_df_copy[~new_df_copy['cluster_accession'].isin(existing_ids)]

    existing_clusters_updated = 0
    if len(new_for_existing_df) > 0:
        meta_indexed = existing_meta.set_index('cluster_id')

        if method_type == 'bin':
            strategy = create_consensus_strategy('bin')
            for cid in new_for_existing_df['cluster_accession'].unique():
                if cid not in meta_indexed.index:
                    continue
                new_spectra = new_for_existing_df[new_for_existing_df['cluster_accession'] == cid]
                existing_row = meta_indexed.loc[cid]
                existing_mz = existing_row['consensus_mz_array']
                existing_int = existing_row['consensus_intensity_array']
                if existing_mz is None or len(existing_mz) == 0:
                    continue

                rows = [{
                    'usi': f'existing_consensus:{cid}',
                    'pepmass': float(existing_row['precursor_mz']),
                    'charge': int(existing_row['charge']),
                    'mz_array': np.array(existing_mz, dtype=np.float32),
                    'intensity_array': np.array(existing_int, dtype=np.float32),
                    'peptidoform': str(existing_row['peptidoform']),
                    'posterior_error_probability': float(existing_row['best_pep']),
                    'global_qvalue': float(existing_row.get('best_qvalue', 1.0)),
                    'cluster_accession': cid,
                }]
                for _, spec in new_spectra.iterrows():
                    pep_str = str(spec['peptidoform'])
                    if '/' not in pep_str:
                        pep_str = pep_str + '/' + str(int(spec['charge']))
                    rows.append({
                        'usi': str(spec['usi']),
                        'pepmass': float(spec['pepmass']),
                        'charge': int(spec['charge']),
                        'mz_array': np.array(spec['mz_array'], dtype=np.float32),
                        'intensity_array': np.array(spec['intensity_array'], dtype=np.float32),
                        'peptidoform': pep_str,
                        'posterior_error_probability': float(spec['posterior_error_probability']),
                        'global_qvalue': float(spec['global_qvalue']),
                        'cluster_accession': cid,
                    })

                mini_df = pd.DataFrame(rows)
                try:
                    consensus_df, _ = strategy.consensus_spectrum_aggregation(mini_df)
                    if not consensus_df.empty:
                        row_out = consensus_df.iloc[0]
                        meta_indexed.loc[cid, 'consensus_mz_array'] = (
                            row_out['mz_array'].astype(np.float32).tolist()
                            if isinstance(row_out['mz_array'], np.ndarray)
                            else list(row_out['mz_array']))
                        meta_indexed.loc[cid, 'consensus_intensity_array'] = (
                            row_out['intensity_array'].astype(np.float32).tolist()
                            if isinstance(row_out['intensity_array'], np.ndarray)
                            else list(row_out['intensity_array']))
                        meta_indexed.loc[cid, 'precursor_mz'] = float(row_out['pepmass'])
                        existing_clusters_updated += 1
                except Exception as e:
                    logger.warning(f"BIN consensus failed for cluster {cid}: {e}")
        else:
            # BEST: swap spectrum when new PSM has better qvalue
            best_new_idx = new_for_existing_df.groupby('cluster_accession')['global_qvalue'].idxmin()
            best_new = new_for_existing_df.loc[best_new_idx].set_index('cluster_accession')
            joined = meta_indexed.loc[meta_indexed.index.isin(best_new.index)].copy()
            joined['new_qvalue'] = best_new['global_qvalue'].astype(float)
            better_ids = joined[joined['new_qvalue'] < joined['best_qvalue']].index

            for cid in better_ids:
                new_row = best_new.loc[cid]
                meta_indexed.loc[cid, 'consensus_mz_array'] = (
                    new_row['mz_array'].astype(np.float32).tolist()
                    if isinstance(new_row['mz_array'], np.ndarray)
                    else list(new_row['mz_array']))
                meta_indexed.loc[cid, 'consensus_intensity_array'] = (
                    new_row['intensity_array'].astype(np.float32).tolist()
                    if isinstance(new_row['intensity_array'], np.ndarray)
                    else list(new_row['intensity_array']))
                meta_indexed.loc[cid, 'precursor_mz'] = float(new_row['pepmass'])
                pep_str = str(new_row['peptidoform'])
                if '/' not in pep_str:
                    pep_str = pep_str + '/' + str(int(new_row['charge']))
                meta_indexed.loc[cid, 'peptidoform'] = pep_str
                meta_indexed.loc[cid, 'peptide_sequence'] = re.sub(
                    r'\[.*?\]', '', pep_str.split('/')[0])
            existing_clusters_updated = len(better_ids)

        existing_meta = meta_indexed.reset_index().rename(columns={'index': 'cluster_id'})
        if 'cluster_id' not in existing_meta.columns:
            existing_meta = existing_meta.reset_index()

    new_clusters_added = 0
    if len(new_brand_new_df) > 0:
        new_meta_df = _build_new_cluster_metadata(
            new_brand_new_df, species, instrument, charge_int, method_type)
        new_clusters_added = len(new_meta_df)
        meta = pd.concat([existing_meta, new_meta_df], ignore_index=True)
    else:
        meta = existing_meta

    meta['cluster_id'] = meta['cluster_id'].astype(str)

    # Update stats from authoritative membership
    stats_map = stats_purity.set_index('cluster_id')
    for col in ['member_count', 'project_count', 'best_pep', 'best_qvalue', 'purity']:
        if col in stats_map.columns:
            meta[col] = meta['cluster_id'].map(stats_map[col]).fillna(meta[col])

    # Provenance columns
    meta['is_reused_cluster'] = meta['cluster_id'].isin(existing_ids)
    sources = _aggregate_source_datasets(out_membership)
    sources_map = dict(zip(sources['cluster_id'].astype(str), sources['source_datasets']))
    meta['source_datasets'] = meta['cluster_id'].map(sources_map).apply(
        lambda x: list(x) if isinstance(x, (list, np.ndarray)) else [project_accession])

    # Ensure correct dtypes
    meta['charge'] = meta['charge'].astype(int)
    meta['member_count'] = meta['member_count'].astype(int)
    meta['project_count'] = meta['project_count'].astype(int)
    meta['purity'] = meta['purity'].astype(float)

    table = pa.Table.from_pandas(
        meta[[f.name for f in CLUSTER_META_SCHEMA]],
        schema=CLUSTER_META_SCHEMA, preserve_index=False)
    pq.write_table(table, out_metadata, compression='zstd')

    logger.info(
        f"Incremental merge complete: {new_clusters_added} new, "
        f"{existing_clusters_updated} updated, {len(meta)} total clusters")
    return out_metadata, out_membership


def _build_new_cluster_metadata(
    new_brand_new_df: pd.DataFrame,
    species: str,
    instrument: str,
    charge_int: int,
    method_type: str,
) -> pd.DataFrame:
    """Build metadata rows for brand-new clusters (not present in existing DB)."""
    from pyspectrafuse.commands.spectrum2msp import create_consensus_strategy

    new_cluster_ids = new_brand_new_df['cluster_accession'].unique()

    if method_type == 'bin' and len(new_brand_new_df) > len(new_cluster_ids):
        strategy = create_consensus_strategy('bin')
        counts = new_brand_new_df['cluster_accession'].value_counts()
        multi_mask = new_brand_new_df['cluster_accession'].isin(counts[counts > 1].index)
        multi_df = new_brand_new_df[multi_mask].copy()
        single_new_df = new_brand_new_df[~multi_mask]

        parts = []
        if not multi_df.empty:
            pep_s = multi_df['peptidoform'].astype(str)
            needs_ch = ~pep_s.str.contains('/', regex=False)
            multi_df['peptidoform'] = pep_s.where(
                ~needs_ch, pep_s + '/' + multi_df['charge'].astype(int).astype(str))
            try:
                consensus_df, singles_df = strategy.consensus_spectrum_aggregation(multi_df)
                for source_df in [consensus_df, singles_df]:
                    if source_df is None or source_df.empty:
                        continue
                    parts.append(_metadata_part(
                        source_df, species, instrument, charge_int, method_type))
            except Exception as e:
                logger.warning(f"BIN consensus failed for new clusters: {e}")
                single_new_df = pd.concat([single_new_df, multi_df])

        if not single_new_df.empty:
            parts.append(_best_per_cluster_metadata_part(
                single_new_df, species, instrument, charge_int, method_type))

        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    return _best_per_cluster_metadata_part(
        new_brand_new_df, species, instrument, charge_int, method_type)


def _best_per_cluster_metadata_part(
    df: pd.DataFrame, species: str, instrument: str,
    charge_int: int, method_type: str,
) -> pd.DataFrame:
    """Pick the best-qvalue spectrum per cluster and build metadata rows."""
    best_idx = df.groupby('cluster_accession')['global_qvalue'].idxmin()
    best = df.loc[best_idx].set_index('cluster_accession')
    pep_series = best['peptidoform'].astype(str)
    needs_charge = ~pep_series.str.contains('/', regex=False)
    pep_series = pep_series.where(
        ~needs_charge,
        pep_series + '/' + best['charge'].astype(int).astype(str))

    return pd.DataFrame({
        'cluster_id': best.index.astype(str),
        'species': species,
        'instrument': instrument,
        'charge': charge_int,
        'peptidoform': pep_series.values,
        'peptide_sequence': pep_series.str.split('/').str[0].str.replace(
            r'\[.*?\]', '', regex=True).values,
        'consensus_mz_array': best['mz_array'].values,
        'consensus_intensity_array': best['intensity_array'].values,
        'consensus_method': method_type,
        'precursor_mz': best['pepmass'].astype(float).values,
        'member_count': 1,
        'project_count': 1,
        'best_pep': best['posterior_error_probability'].astype(float).values,
        'best_qvalue': best['global_qvalue'].astype(float).values,
        'purity': 1.0,
    })


def _metadata_part(
    source_df: pd.DataFrame, species: str, instrument: str,
    charge_int: int, method_type: str,
) -> pd.DataFrame:
    """Materialize a metadata block from a consensus-strategy output df."""
    cid_source = (source_df['cluster_accession'].astype(str).values
                  if 'cluster_accession' in source_df.columns
                  else source_df.index.astype(str).values)
    return pd.DataFrame({
        'cluster_id': cid_source,
        'species': species,
        'instrument': instrument,
        'charge': charge_int,
        'peptidoform': source_df['peptidoform'].astype(str).values,
        'peptide_sequence': source_df['peptidoform'].astype(str).str.split('/').str[0].str.replace(
            r'\[.*?\]', '', regex=True).values,
        'consensus_mz_array': source_df['mz_array'].values,
        'consensus_intensity_array': source_df['intensity_array'].values,
        'consensus_method': method_type,
        'precursor_mz': pd.to_numeric(source_df['pepmass'], errors='coerce').astype(float).values,
        'member_count': 1,
        'project_count': 1,
        'best_pep': pd.to_numeric(
            source_df['posterior_error_probability'], errors='coerce').astype(float).values,
        'best_qvalue': pd.to_numeric(
            source_df.get('global_qvalue', source_df['posterior_error_probability']),
            errors='coerce').astype(float).values,
        'purity': 1.0,
    })
