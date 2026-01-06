from typing import Tuple
from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
from pyspectrafuse.common.msp_utils import MspUtil


class BestSpectrumStrategy(ConsensusStrategy):
    """Best spectrum strategy selects the spectrum with the best metric value.
    
    For clusters with multiple spectra, selects the spectrum with the minimum
    value of the specified filter metric (typically posterior_error_probability
    or global_qvalue) as the consensus spectrum.
    """

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, 
                                       filter_metrics: str = 'global_qvalue') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Generate consensus spectrum by selecting best spectrum per cluster.
        
        Args:
            cluster_df: DataFrame containing cluster spectrum data
            filter_metrics: Metric column name to use for selection
            
        Returns:
            Tuple of (consensus_spectrum_df, single_spectrum_df)
        """
        return self.classify_cluster_group(cluster_df, filter_metrics)

    def classify_cluster_group(self, df: pd.DataFrame, filter_metrics: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Classify clusters and select best spectrum for each multi-spectrum cluster.
        
        Args:
            df: DataFrame containing cluster data
            filter_metrics: Metric column name for selecting best spectrum
            
        Returns:
            Tuple of (best_spectrum_df, single_spectrum_df)
        """
        counts = self.get_cluster_counts(df)
        counts_dict = counts.to_dict()
        
        # Use vectorized map instead of apply for better performance
        df['Nreps'] = df['cluster_accession'].map(counts_dict)
        df['peptidoform'] = df['peptidoform'] + '/' + df['charge'].astype(str)
        
        # For non-single clusters, select the spectrum with the smallest metric value as consensus
        count_greater_than_1 = df[np.isin(df['cluster_accession'], counts[counts > 1].index)]
        count_greater_than_1_groups = count_greater_than_1.groupby('cluster_accession')
        
        # Select spectrum with minimum metric value
        best_spectrum_group = count_greater_than_1_groups.apply(
            self.top_n_rows, column=filter_metrics, n=1, include_groups=False)
        best_spectrum_df = best_spectrum_group.reset_index(drop=True)

        # Single spectrum clusters
        single_spectrum_df = df[np.isin(df['cluster_accession'], counts[counts == 1].index)]

        return best_spectrum_df, single_spectrum_df

    @staticmethod
    def _create_msp_dict_row(row: pd.Series, msp_util: MspUtil, nreps: int) -> dict:
        """Create MSP dictionary for a single row.
        
        Args:
            row: Pandas Series with spectrum data
            msp_util: MspUtil instance
            nreps: Number of replicates
            
        Returns:
            MSP format dictionary
        """
        return MspUtil.get_msp_dict(
            name=row['usi'].split(':')[-1],
            mw=row['pepmass'],
            num_peaks=row['mz_array'].shape[0],
            comment=f'clusterID={msp_util.usi_to_uuid([row["usi"]])} Nreps={nreps} PEP={row["posterior_error_probability"]}',
            mz_arr=row['mz_array'],
            intensity_arr=row['intensity_array']
        )

    @staticmethod
    def single_and_consensus_to_msp_dict(single_df: pd.DataFrame, 
                                         consensus_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Convert consensus and single spectra DataFrames to MSP dictionaries.
        
        Args:
            single_df: DataFrame with single spectra
            consensus_df: DataFrame with consensus spectra
            
        Returns:
            Tuple of (consensus_df, single_df) with added 'Msp_dict' column
        """
        msp_util = MspUtil()
        
        # Use helper method to reduce duplication
        consensus_df['Msp_dict'] = consensus_df.apply(
            lambda row: BestSpectrumStrategy._create_msp_dict_row(row, msp_util, row['Nreps']), axis=1)
        
        single_df['Msp_dict'] = single_df.apply(
            lambda row: BestSpectrumStrategy._create_msp_dict_row(row, msp_util, 1), axis=1)

        return consensus_df, single_df




