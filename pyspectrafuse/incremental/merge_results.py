"""Merge incremental clustering results into existing parquet tables.

After resolve_clusters maps new cluster IDs to existing/fresh ones, this
module appends new PSMs to psm_cluster_membership.parquet and regenerates
cluster_metadata.parquet with updated statistics and consensus spectra.

Uses DuckDB for aggregation and I/O operations.
"""
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pyspectrafuse.common.schemas import CLUSTER_META_SCHEMA, PSM_MEMBERSHIP_SCHEMA
from pyspectrafuse.common.duckdb_ops import (
    append_membership_dedup,
    compute_stats_and_purity,
    write_parquet_from_df,
)

logger = logging.getLogger(__name__)


def append_psm_membership(
    existing_path: str,
    new_psms_df: pd.DataFrame,
    output_path: Optional[str] = None,
) -> str:
    """Append new PSM rows to the membership parquet, deduplicating by USI.

    Args:
        existing_path: Path to existing psm_cluster_membership.parquet.
        new_psms_df: DataFrame with same schema as the membership table.
            Must have 'resolved_cluster_id' which will be renamed to 'cluster_id'.
        output_path: Where to write. Defaults to overwriting existing_path.

    Returns:
        Output path.
    """
    if output_path is None:
        output_path = existing_path

    # Normalize column name
    if 'resolved_cluster_id' in new_psms_df.columns:
        new_psms_df = new_psms_df.rename(columns={'resolved_cluster_id': 'cluster_id'})

    return append_membership_dedup(existing_path, new_psms_df, output_path)


def rebuild_cluster_metadata(
    psm_membership_path: str,
    consensus_df: pd.DataFrame,
    single_df: pd.DataFrame,
    species: str,
    instrument: str,
    charge: str,
    method_type: str,
    output_path: str,
) -> str:
    """Rebuild cluster_metadata.parquet from the full PSM membership table.

    This recomputes member_count, project_count, purity, best_pep, and
    best_qvalue from the authoritative PSM table, then joins with consensus
    spectra.

    Args:
        psm_membership_path: Path to the (already merged) PSM membership parquet.
        consensus_df: Consensus spectra DataFrame (from consensus strategy).
        single_df: Single-spectrum clusters DataFrame.
        species, instrument, charge, method_type: Metadata fields.
        output_path: Where to write the new cluster_metadata.parquet.

    Returns:
        Output path.
    """
    # Compute stats and purity via DuckDB (reads parquet directly, no pandas load)
    stats_purity = compute_stats_and_purity(psm_membership_path)

    # Build metadata from consensus + single spectra (stays in pandas for list columns)
    charge_int = int(str(charge).replace('charge', ''))
    parts = []

    for source_df in [consensus_df, single_df]:
        if source_df is None or source_df.empty:
            continue

        part = pd.DataFrame()
        part['cluster_id'] = source_df['cluster_accession'].astype(str).values
        part['species'] = species
        part['instrument'] = instrument
        part['charge'] = charge_int
        part['peptidoform'] = source_df['peptidoform'].astype(str).values
        part['peptide_sequence'] = (
            part['peptidoform'].str.split('/').str[0]
            .str.replace(r'\[.*?\]', '', regex=True))

        def _to_f32_list(arr):
            if isinstance(arr, np.ndarray):
                return arr.astype(np.float32).tolist()
            if isinstance(arr, (list, tuple)):
                return [float(x) for x in arr]
            return []

        part['consensus_mz_array'] = source_df['mz_array'].apply(_to_f32_list).values
        part['consensus_intensity_array'] = source_df['intensity_array'].apply(_to_f32_list).values
        part['consensus_method'] = method_type
        part['precursor_mz'] = pd.to_numeric(
            source_df['pepmass'], errors='coerce').astype(float).values
        parts.append(part)

    if not parts:
        logger.warning("No consensus spectra to write")
        return output_path

    meta_df = pd.concat(parts, ignore_index=True)

    # Merge with stats and purity from DuckDB
    meta_df = meta_df.merge(stats_purity, on='cluster_id', how='left')
    meta_df['member_count'] = meta_df['member_count'].fillna(1).astype(int)
    meta_df['project_count'] = meta_df['project_count'].fillna(1).astype(int)
    meta_df['purity'] = meta_df['purity'].fillna(1.0)

    # Write with canonical schema
    table = pa.Table.from_pandas(meta_df, schema=CLUSTER_META_SCHEMA, preserve_index=False)
    pq.write_table(table, output_path, compression='zstd')

    logger.info(f"Cluster metadata: {len(meta_df)} clusters written to {output_path}")
    return output_path
