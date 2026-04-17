"""Centralized Parquet schemas for cluster DB output files.

These schemas define the canonical column names and types for the two output
parquet files produced by the clustering pipeline:
- cluster_metadata.parquet: one row per cluster (consensus spectrum + stats)
- psm_cluster_membership.parquet: one row per PSM (cluster assignment + scores)
"""
import pyarrow as pa

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
