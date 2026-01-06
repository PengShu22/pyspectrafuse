from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
from pyspectrafuse.common.msp_utils import MspUtil


class BestSpectrumStrategy(ConsensusStrategy):
    """
    only consider the posterior_error_probability value
    """

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, filter_metrics='global_qvalue'):
        return self.classify_cluster_group(cluster_df, filter_metrics)

    def classify_cluster_group(self, df: pd.DataFrame, filter_metrics):
        counts = self.get_cluster_counts(df)
        counts_dict = counts.to_dict()
        df['Nreps'] = df['cluster_accession'].apply(lambda x: counts_dict[x])
        df['peptidoform'] = df['peptidoform'] + '/' + str(df['charge'].to_list()[0])
        # For non-single clusters, select the spectrum with the smallest posterior_error_probability as consensus
        count_greater_than_2 = df[np.in1d(df['cluster_accession'], counts[counts > 1].index)]
        count_greater_than_2_groups = count_greater_than_2.groupby('cluster_accession')
        # Select spectrum with minimum posterior_error_probability
        best_spectrum_group = count_greater_than_2_groups.apply(self.top_n_rows, column=filter_metrics,
                                                                n=1)

        best_spectrum_df = best_spectrum_group.reset_index(drop=True)  # Generated consensus spectrum

        single_spectrum_df = df[
            np.in1d(df['cluster_accession'], counts[counts == 1].index)]  # spectrum number equal 1

        return best_spectrum_df, single_spectrum_df

    @staticmethod
    def single_and_consensus_to_msp_dict(single_df: pd.DataFrame, consensus_df: pd.DataFrame):
        mspUtil = MspUtil()
        consensus_df['Msp_dict'] = consensus_df.apply(lambda row: MspUtil.get_msp_dict(
            name=row['usi'].split(':')[-1],
            mw=row['pepmass'],
            num_peaks=row['mz_array'].shape[0],
            comment=f'clusterID={mspUtil.usi_to_uuid([row["usi"]])} Nreps={row["Nreps"]} PEP={row["posterior_error_probability"]}',
            mz_arr=row['mz_array'],
            intensity_arr=row['intensity_array']
        ), axis=1)

        single_df['Msp_dict'] = single_df.apply(lambda row: MspUtil.get_msp_dict(
            name=row['usi'].split(':')[-1],
            mw=row['pepmass'],
            num_peaks=row['mz_array'].shape[0],
            comment=f'clusterID={mspUtil.usi_to_uuid([row["usi"]])} Nreps=1 PEP={row["posterior_error_probability"]}',
            mz_arr=row['mz_array'],
            intensity_arr=row['intensity_array']
        ), axis=1)

        return consensus_df, single_df




