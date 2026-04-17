"""Resolve incremental clustering results using scan_titles (dat-mode).

Replaces resolve_clusters.py which parsed MGF files. This version reads
scan_titles.txt files to identify representative spectra (titles starting
with ``rep:``) and map new MaRaCluster cluster IDs back to existing ones.

Terminology:
- *old cluster*: a cluster from the previous round (has a ``cluster_id`` in
  ``cluster_metadata.parquet``).
- *representative*: the consensus spectrum written to .dat with scan_title
  ``rep:{cluster_id}``.
- *new cluster*: a cluster assigned by MaRaCluster in the incremental round.
- *resolved cluster*: the final cluster ID — reuses the old ID when a
  representative is present, or mints a new UUID when no representative maps.
"""
import logging
import re
import uuid
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

_REP_PREFIX = "rep:"


def parse_cluster_tsv(tsv_path: str) -> pd.DataFrame:
    """Read a MaRaCluster TSV into a DataFrame.

    Returns:
        DataFrame with columns [mgf_path, scannr, new_cluster_id].
    """
    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=['mgf_path', 'scannr', 'new_cluster_id'])
    df.dropna(inplace=True)
    df['scannr'] = df['scannr'].astype(int)
    # Strip window suffix from MGF basename for matching with scan_titles
    df['mgf_basename'] = df['mgf_path'].apply(
        lambda p: re.sub(r'_w\d+\.mgf$', '.mgf', Path(p).name))
    return df


def parse_all_scan_titles(scan_titles_files: List[str]) -> pd.DataFrame:
    """Parse scan_titles.txt files into a unified lookup table.

    Returns:
        DataFrame with columns [mgf_basename, scannr, title, is_rep, old_cluster_id].
    """
    rows = []
    for titles_path in scan_titles_files:
        # Derive mgf_basename from the scan_titles filename
        # e.g., "PXD014877.psm_charge2.scan_titles.txt" -> "PXD014877.psm_charge2.mgf"
        stem = Path(titles_path).name.replace('.scan_titles.txt', '')
        mgf_basename = f'{stem}.mgf'

        with open(titles_path) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t', 2)
                if len(parts) < 3:
                    continue
                scannr = int(parts[1])
                title = parts[2]
                is_rep = title.startswith(_REP_PREFIX)
                old_cluster_id = title[len(_REP_PREFIX):] if is_rep else None
                rows.append((mgf_basename, scannr, title, is_rep, old_cluster_id))

    return pd.DataFrame(rows, columns=[
        'mgf_basename', 'scannr', 'title', 'is_rep', 'old_cluster_id'])


def resolve_incremental_clusters(
    cluster_tsv_path: str,
    scan_titles_files: List[str],
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Map new MaRaCluster cluster IDs to resolved (old or fresh) IDs.

    Algorithm:
    1. Parse the MaRaCluster TSV to get (mgf_basename, scannr) → new_cluster_id.
    2. Join with scan_titles to identify which spectra are representatives.
    3. For representatives, map new_cluster_id → old_cluster_id.
    4. For new clusters with no representative, mint a fresh UUID.

    If multiple old clusters merge into one new cluster (representatives from
    different old clusters land in the same new cluster), the first old cluster
    ID wins (deterministic).

    Args:
        cluster_tsv_path: Path to MaRaCluster output TSV.
        scan_titles_files: Paths to all scan_titles.txt files (reps + new data).

    Returns:
        Tuple of:
        - id_map: {new_cluster_id → resolved_cluster_id}
        - new_spectra_df: DataFrame of non-representative spectra with
          columns [mgf_basename, scannr, new_cluster_id, resolved_cluster_id, title].
    """
    cluster_df = parse_cluster_tsv(cluster_tsv_path)
    titles_df = parse_all_scan_titles(scan_titles_files)

    logger.info(f"Cluster TSV: {len(cluster_df)} rows")
    logger.info(f"Scan titles: {len(titles_df)} entries "
                f"({titles_df['is_rep'].sum()} representatives)")

    # Join cluster assignments with scan_titles
    merged = cluster_df.merge(titles_df, on=['mgf_basename', 'scannr'], how='left')

    # Split into representative and new spectra
    rep_rows = merged[merged['is_rep'] == True].copy()
    new_rows = merged[merged['is_rep'] != True].copy()

    logger.info(f"Matched: {len(rep_rows)} representatives, {len(new_rows)} new spectra")

    # Build new_cluster_id → old_cluster_id mapping from representatives
    new_to_old: Dict[str, list] = {}
    for _, row in rep_rows.iterrows():
        new_cid = str(row['new_cluster_id'])
        old_cid = row['old_cluster_id']
        if old_cid is not None:
            new_to_old.setdefault(new_cid, []).append(old_cid)

    # Resolve: for each new_cluster_id, pick the old_id (or mint fresh)
    id_map: Dict[str, str] = {}
    merge_count = 0

    for new_cid in cluster_df['new_cluster_id'].unique():
        new_cid_str = str(new_cid)
        old_ids = new_to_old.get(new_cid_str, [])
        if len(old_ids) == 0:
            id_map[new_cid_str] = str(uuid.uuid4())
        elif len(old_ids) == 1:
            id_map[new_cid_str] = old_ids[0]
        else:
            id_map[new_cid_str] = old_ids[0]
            merge_count += 1
            logger.info(f"Cluster merge: new {new_cid_str} <- old {old_ids}")

    if merge_count > 0:
        logger.info(f"Total cluster merges: {merge_count}")

    # Tag new spectra with resolved IDs
    new_rows = new_rows.copy()
    new_rows['resolved_cluster_id'] = (
        new_rows['new_cluster_id'].astype(str).map(id_map))

    return id_map, new_rows
