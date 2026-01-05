import logging
import numpy as np
import pandas as pd
import spectrum_utils.spectrum as sus

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


# 定义共识谱生成策略抽象基类
class ConsensusStrategy:
    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, filter_metrics):
        pass

    @staticmethod
    def top_n_rows(df, column, n):
        # 获取列的值
        values = df[column].values
        # 计算前n个最小值的索引
        indices = np.argpartition(values, n)[:n]
        # 返回对应的行
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
    def get_cluster_counts(df: pd.DataFrame):
        return df['cluster_accession'].value_counts(sort=False)

    def classify_cluster_group(self, df: pd.DataFrame, filter_metrics):
        counts = self.get_cluster_counts(df)
        counts_dict = counts.to_dict()
        df['Nreps'] = df['cluster_accession'].apply(lambda x: counts_dict[x])
        df['peptidoform'] = df['peptidoform'] + '/' + str(df['charge'].to_list()[0])
        count_greater_than_10 = df[
            np.in1d(df['cluster_accession'], counts[counts > 10].index)]
        count_greater_than_10_groups = count_greater_than_10.groupby('cluster_accession')

        count_top10 = count_greater_than_10_groups.apply(self.top_n_rows, column=filter_metrics, n=10)  # top 10

        count_median_10 = df[np.in1d(df['cluster_accession'],
                                     counts[(counts >= 2) & (
                                             counts <= 10)].index)]  # spectrum between 2 and 10
        # 合并过滤后的top10和median
        merge_df = pd.DataFrame(np.vstack([count_top10.values, count_median_10.values]),
                                columns=count_top10.columns)

        single_df = df[
            np.in1d(df['cluster_accession'], counts[counts == 1].index)]  # spectrum number equal 1

        return merge_df, single_df

