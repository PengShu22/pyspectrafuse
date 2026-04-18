from enum import Enum


class UseCol(Enum):
    """Column names used for parquet file operations.

    Supports both legacy msnet format and QPX format by reading all possible
    column names and normalizing at load time via ParquetSchemaAdapter.
    """
    # Columns needed for MGF conversion (legacy msnet format)
    PARQUET_COL_TO_MGF = ['sequence', 'mz_array', 'intensity_array',
                          'precursor_charge', 'exp_mass_to_charge', 'scan', 'reference_file_name']

    # Columns needed for MSP/cluster-parquet conversion
    # We list both legacy and QPX names; the adapter picks whichever exists
    PARQUET_COL_TO_MSP = ['sequence', 'peptidoform',
                          'mz_array', 'intensity_array',
                          'precursor_charge', 'charge',  # legacy / QPX
                          'exp_mass_to_charge', 'observed_mz',  # legacy / QPX
                          'reference_file_name', 'run_file_name',  # legacy / QPX
                          'scan']

    PARQUET_COL_TO_FILTER = ['posterior_error_probability', 'global_qvalue']

    PARQUET_COL_TO_USI = ['reference_file_name', 'run_file_name',
                          'scan', 'sequence', 'precursor_charge', 'charge']


# Canonical column names used internally by pyspectrafuse
# All input formats are normalized to these names
CANONICAL_CHARGE = 'charge'
CANONICAL_OBSERVED_MZ = 'pepmass'
CANONICAL_RUN_FILE = 'reference_file_name'
CANONICAL_SCAN = 'scan'


class ParquetSchemaAdapter:
    """Adapts parquet DataFrames from legacy msnet or QPX format to canonical names.

    Call adapt(df) after reading a parquet batch to normalize column names.
    """

    # (canonical_name, [candidate_names_in_priority_order])
    COLUMN_MAP = [
        ('charge', ['charge', 'precursor_charge']),
        ('pepmass', ['observed_mz', 'exp_mass_to_charge']),
        ('reference_file_name', ['run_file_name', 'reference_file_name']),
    ]

    @staticmethod
    def available_columns(parquet_columns: list, needed: list) -> list:
        """Return only the columns from `needed` that exist in the parquet file."""
        parquet_set = set(parquet_columns)
        return [c for c in needed if c in parquet_set]

    @classmethod
    def adapt(cls, df):
        """Rename columns to canonical names. Handles both msnet and QPX schemas."""
        for canonical, candidates in cls.COLUMN_MAP:
            if canonical in df.columns:
                continue
            for candidate in candidates:
                if candidate in df.columns:
                    df = df.rename(columns={candidate: canonical})
                    break

        # QPX scan is list<int32>, legacy is int32. Normalize to int.
        if 'scan' in df.columns and len(df) > 0:
            import numpy as np
            sample = df['scan'].iloc[0]
            if isinstance(sample, (list, tuple, np.ndarray)):
                df['scan'] = df['scan'].apply(
                    lambda x: int(x[0]) if hasattr(x, '__len__') and len(x) > 0 else 0)

        # If global_qvalue is missing (QPX stores it in additional_scores), extract it
        if 'global_qvalue' not in df.columns and 'additional_scores' in df.columns:
            df['global_qvalue'] = df['additional_scores'].apply(cls._extract_global_qvalue)

        return df

    @staticmethod
    def _extract_global_qvalue(scores):
        """Extract global_qvalue from QPX additional_scores list/array."""
        if scores is None:
            return None
        try:
            for score in scores:
                if isinstance(score, dict) and score.get('score_name') == 'global_qvalue':
                    return score.get('score_value')
        except (TypeError, ValueError):
            pass
        return None


# Magic numbers and constants
class ClusterConstants:
    """Constants used for cluster classification."""
    # Cluster size thresholds
    LARGE_CLUSTER_THRESHOLD = 10  # Clusters with more than 10 spectra
    MEDIUM_CLUSTER_MIN = 2  # Minimum for medium clusters
    MEDIUM_CLUSTER_MAX = 10  # Maximum for medium clusters
    SINGLE_CLUSTER_SIZE = 1  # Single spectrum clusters

    # Top N selection
    TOP_N_FOR_LARGE_CLUSTERS = 10  # Select top 10 spectra from large clusters
    TOP_N_FOR_CONSENSUS = 1  # Select 1 spectrum for consensus
