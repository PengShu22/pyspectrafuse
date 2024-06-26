from typing import Tuple
from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
import logging
import functools
from pyspectrafuse.consensus_strategy.metrics import dot

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


class MostSimilarStrategy(ConsensusStrategy):

    def __init__(self, sim='dot', fragment_mz_tolerance=0.02):
        self.fragment_mz_tolerance = 0.02
        self.sim = sim
        if self.sim == 'dot':
            self.compare_spectra = functools.partial(
                dot, fragment_mz_tolerance=fragment_mz_tolerance)
        else:
            raise ValueError("Unknown spectrum similarity method")

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame,
                                       filter_metrics='global_qvalue') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        该方法对簇中每一个谱图与其簇内其他谱图进行相似距离计算，并对相似距离进行求和，最终选取相似距离之和最大的谱图作为代表谱图。
        对于出现多个距离之和相同的情况，MOST总是随机选择一个谱图作为代表。该算法所生成的一致性谱图来源于谱图本身，
        和BEST算法类似，并不会生成新的谱图。
        """
        if self.sim == 'dot':
            merge_median_and_top, single_spectrum_df = self.classify_cluster_group(cluster_df, filter_metrics)

            merge_median_and_top['ms2spectrumObj'] = merge_median_and_top.apply(
                lambda row: self.get_Ms2SpectrumObj(row),
                axis=1)
            # 对accession进行分组，并对每个组应用print_dict函数
            res = merge_median_and_top.groupby('cluster_accession').apply(self._select_representative)
            logger.info(
                f"consensus spectrum shape: {res.shape} ; single spectrum shape: {single_spectrum_df.shape}")

            return res, single_spectrum_df

    def _select_representative(self, single_group: pd.DataFrame) -> pd.Series:
        # 内部应该为一个簇生成一个字典，key为usi， value为sus.MsmsSpectrum对象
        cluster_spectra = single_group.set_index('usi')['ms2spectrumObj'].to_dict()
        spectra_keys = list(cluster_spectra.keys())

        sim_matrix = np.zeros((len(spectra_keys), len(spectra_keys)))
        for i in range(len(spectra_keys)):
            for j in range(i, len(spectra_keys)):
                # 传入两个谱图 使用dot计算相似程度 两个谱图可以相同
                sim_matrix[i, j] = sim_matrix[j, i] = self.compare_spectra(
                    cluster_spectra[spectra_keys[i]],
                    cluster_spectra[spectra_keys[j]])
        # Find the spectrum with the maximum similarity to all other spectra.
        max_sim_index = sim_matrix.sum(axis=0).argmax()
        max_sim_spectrum = cluster_spectra[spectra_keys[max_sim_index]]
        max_sim_usi = spectra_keys[max_sim_index]

        # 找到最相似谱图对应的行并返回
        max_sim_row = single_group[single_group['usi'] == max_sim_usi]
        return max_sim_row
