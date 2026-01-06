import logging
from typing import Tuple
import numpy as np
import pandas as pd
import spectrum_utils.spectrum as sus
from pyspectrafuse.common.constant import ClusterConstants

logger = logging.getLogger(__name__)


# Abstract base class for consensus spectrum generation strategies
class ConsensusStrategy:
    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, 
                                       filter_metrics: str = 'global_qvalue') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Generate consensus spectrum from cluster dataframe.
        
        Args:
            cluster_df: DataFrame containing cluster spectrum data
            filter_metrics: Metric to use for filtering (default: 'global_qvalue')
            
        Returns:
            Tuple of (consensus_spectrum_df, single_spectrum_df)
        """
        raise NotImplementedError("Subclasses must implement consensus_spectrum_aggregation")

    @staticmethod
    def top_n_rows(df: pd.DataFrame, column: str, n: int) -> pd.DataFrame:
        """Get the top n rows with smallest values in the specified column.
        
        Args:
            df: DataFrame to filter
            column: Column name to sort by
            n: Number of rows to return
            
        Returns:
            DataFrame with top n rows
        """
        # Get column values
        values = df[column].values
        # Calculate indices of the n smallest values
        indices = np.argpartition(values, n)[:n]
        # Return corresponding rows
        return df.iloc[indices]

    @staticmethod
    def get_Ms2SpectrumObj(row: pd.Series) -> sus.MsmsSpectrum:
        spectrum = sus.MsmsSpectrum(
            row['usi'],
            row['pepmass'],
            row['charge'],
            row['mz_array'],
            row['intensity_array']
        )

        return spectrum

    @staticmethod
    def get_Ms2SpectrumDict(row: pd.Series) -> dict:
        res_dict = {
            "params": {
                'pepmass': [row['pepmass']],
                'charge': [row['charge']]
            },
            'm/z array': row['mz_array'],
            'intensity array': row['intensity_array']
        }

        return res_dict

    @staticmethod
    def get_cluster_counts(df: pd.DataFrame) -> pd.Series:
        """Get count of spectra per cluster.
        
        Args:
            df: DataFrame with cluster_accession column
            
        Returns:
            Series with cluster counts
        """
        return df['cluster_accession'].value_counts(sort=False)

    def classify_cluster_group(self, df: pd.DataFrame, filter_metrics: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Classify clusters into groups based on size and filter metrics.
        
        Args:
            df: DataFrame containing cluster data
            filter_metrics: Metric column name for filtering
            
        Returns:
            Tuple of (merged_df, single_df) for clusters with >1 and ==1 spectra
        """
        counts = self.get_cluster_counts(df)
        counts_dict = counts.to_dict()
        
        # Use vectorized map instead of apply for better performance
        df['Nreps'] = df['cluster_accession'].map(counts_dict)
        df['peptidoform'] = df['peptidoform'] + '/' + df['charge'].astype(str)
        
        # Large clusters (>10 spectra): select top N
        large_clusters = counts[counts > ClusterConstants.LARGE_CLUSTER_THRESHOLD]
        count_greater_than_10 = df[np.isin(df['cluster_accession'], large_clusters.index)]
        count_greater_than_10_groups = count_greater_than_10.groupby('cluster_accession')
        count_top10 = count_greater_than_10_groups.apply(
            self.top_n_rows, column=filter_metrics, 
            n=ClusterConstants.TOP_N_FOR_LARGE_CLUSTERS, include_groups=False)

        # Medium clusters (2-10 spectra): use all
        medium_clusters = counts[(counts >= ClusterConstants.MEDIUM_CLUSTER_MIN) & 
                                (counts <= ClusterConstants.MEDIUM_CLUSTER_MAX)]
        count_median_10 = df[np.isin(df['cluster_accession'], medium_clusters.index)]
        
        # Merge filtered top10 and median results
        merge_df = pd.DataFrame(np.vstack([count_top10.values, count_median_10.values]),
                                columns=count_top10.columns)

        # Single spectrum clusters
        single_df = df[np.isin(df['cluster_accession'], 
                               counts[counts == ClusterConstants.SINGLE_CLUSTER_SIZE].index)]

        return merge_df, single_df

