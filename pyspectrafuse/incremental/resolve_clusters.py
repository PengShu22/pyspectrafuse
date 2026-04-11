"""Resolve incremental clustering results back to existing cluster IDs.

After MaRaCluster runs on representatives + new spectra, new cluster IDs are
assigned. This module maps those back to the original cluster IDs by tracing
which representative ended up in which new cluster.

Terminology:
- *old cluster*: a cluster from the previous round (has a ``cluster_id`` in
  ``cluster_metadata.parquet``).
- *representative*: the consensus spectrum written to MGF with TITLE=``rep:{cluster_id}``.
- *new cluster*: a cluster assigned by MaRaCluster in the incremental round.
- *resolved cluster*: the final cluster ID — reuses the old ID when a
  representative is present, or mints a new UUID when no representative maps.
"""
import logging
import uuid
from typing import Dict, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# Prefix used in representative MGF TITLE fields
_REP_PREFIX = "rep:"


def parse_cluster_tsv(tsv_path: str) -> pd.DataFrame:
    """Read a MaRaCluster TSV into a DataFrame.

    Returns:
        DataFrame with columns [mgf_path, spectrum_index, new_cluster_id].
    """
    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=['mgf_path', 'spectrum_index', 'new_cluster_id'])
    df.dropna(inplace=True)
    df['spectrum_index'] = df['spectrum_index'].astype(int)
    return df


def build_representative_index(rep_mgf_path: str) -> Dict[Tuple[str, int], str]:
    """Scan a representative MGF file and build (mgf_path, index) → old_cluster_id.

    This is lightweight: we only parse TITLE lines, not peak lists.

    Returns:
        Mapping from (mgf_file_basename, spectrum_order) to old cluster_id.
    """
    index: Dict[Tuple[str, int], str] = {}
    current_order = -1  # will be incremented to 0 on first BEGIN IONS

    with open(rep_mgf_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                current_order += 1
            elif line.startswith("TITLE=") and line[6:].startswith(_REP_PREFIX):
                old_cluster_id = line[6 + len(_REP_PREFIX):]
                # MaRaCluster identifies spectra by (file, index)
                # The file is the MGF path itself; index is the 0-based order
                index[current_order] = old_cluster_id

    logger.info(f"Built representative index with {len(index)} entries from {rep_mgf_path}")
    return index


def resolve_incremental_clusters(
    cluster_tsv_path: str,
    rep_mgf_path: str,
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Map new MaRaCluster cluster IDs to resolved (old or fresh) IDs.

    Algorithm:
    1. Parse the MaRaCluster TSV to get (mgf_path, index) → new_cluster_id.
    2. For each representative spectrum, look up which new_cluster_id it was
       assigned to. This gives new_cluster_id → old_cluster_id.
    3. For new clusters that contain NO representative, mint a fresh UUID.
    4. Build a mapping: new_cluster_id → resolved_cluster_id.

    If multiple old clusters merge into one new cluster (representatives from
    different old clusters land in the same new cluster), the old cluster with
    the most members wins; the others are treated as merges.

    Args:
        cluster_tsv_path: Path to MaRaCluster output TSV.
        rep_mgf_path: Path to the representative MGF used in clustering.

    Returns:
        Tuple of:
        - id_map: {new_cluster_id → resolved_cluster_id}
        - new_spectra_df: DataFrame of non-representative spectra with
          columns [mgf_path, spectrum_index, new_cluster_id, resolved_cluster_id].
    """
    cluster_df = parse_cluster_tsv(cluster_tsv_path)
    rep_index = build_representative_index(rep_mgf_path)

    # Identify which rows in the cluster TSV are representatives
    # Representatives came from rep_mgf_path, so their mgf_path will contain
    # the representative MGF filename
    rep_mgf_basename = rep_mgf_path.rsplit('/', 1)[-1] if '/' in rep_mgf_path else rep_mgf_path

    rep_mask = cluster_df['mgf_path'].str.contains(rep_mgf_basename, regex=False)
    rep_rows = cluster_df[rep_mask].copy()
    new_rows = cluster_df[~rep_mask].copy()

    logger.info(f"Cluster TSV: {len(cluster_df)} rows, "
                f"{len(rep_rows)} representatives, {len(new_rows)} new spectra")

    # Build new_cluster_id → old_cluster_id mapping from representatives
    # A new cluster may contain multiple representatives (cluster merge)
    new_to_old: Dict[str, list] = {}
    for _, row in rep_rows.iterrows():
        order = int(row['spectrum_index'])
        old_id = rep_index.get(order)
        if old_id is None:
            logger.warning(f"Representative at index {order} not found in rep index")
            continue
        new_cid = str(row['new_cluster_id'])
        new_to_old.setdefault(new_cid, []).append(old_id)

    # Resolve: for each new_cluster_id, pick the old_id (or mint fresh)
    id_map: Dict[str, str] = {}
    merge_count = 0

    for new_cid in cluster_df['new_cluster_id'].unique():
        new_cid_str = str(new_cid)
        old_ids = new_to_old.get(new_cid_str, [])
        if len(old_ids) == 0:
            # Entirely new cluster — mint a UUID
            id_map[new_cid_str] = str(uuid.uuid4())
        elif len(old_ids) == 1:
            # Clean 1:1 mapping
            id_map[new_cid_str] = old_ids[0]
        else:
            # Merge: multiple old clusters collapsed. Keep first (arbitrary but
            # deterministic); the Nextflow workflow can log merge events.
            id_map[new_cid_str] = old_ids[0]
            merge_count += 1
            logger.info(f"Cluster merge: new {new_cid_str} ← old {old_ids}")

    if merge_count > 0:
        logger.info(f"Total cluster merges: {merge_count}")

    # Tag new (non-representative) spectra with resolved IDs
    new_rows['resolved_cluster_id'] = (
        new_rows['new_cluster_id'].astype(str).map(id_map))

    return id_map, new_rows
