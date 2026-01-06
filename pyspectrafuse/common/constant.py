from enum import Enum


class UseCol(Enum):
    """Column names used for parquet file operations."""
    PARQUET_COL_TO_MGF = ['USI', 'sequence',
                          'mz_array', 'intensity_array', 'charge', 'exp_mass_to_charge']

    PARQUET_COL_TO_MSP = ['USI', 'sequence', 'peptidoform',
                          'mz_array', 'intensity_array', 'charge', 'exp_mass_to_charge']

    PARQUET_COL_TO_FILTER = ['posterior_error_probability', 'global_qvalue']

    PARQUET_COL_TO_USI = ['reference_file_name', 'scan_number', 'sequence', 'charge']


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

