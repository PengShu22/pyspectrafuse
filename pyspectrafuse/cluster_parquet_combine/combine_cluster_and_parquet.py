import logging
import pandas as pd
import numpy as np
from pathlib import Path
import pyarrow.parquet as pq
from collections import defaultdict
from pyspectrafuse.common.constant import UseCol
from pyspectrafuse.common.sdrf_utils import SdrfUtil
from pyspectrafuse.mgf_convert.parquet2mgf import Parquet2Mgf
from pyspectrafuse.common.parquet_utils import ParquetPathHandler

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


class CombineCluster2Parquet:

    def __init__(self):
        self.combine_info_col = ['mz_array', 'intensity_array', 'charge', 'peptidoform',
                                 'pepmass', 'posterior_error_probability', 'global_qvalue']

    @staticmethod
    def map_strategy_to_cluster(map_dict: dict, df_1: pd.DataFrame, classify_path: str, order_range: range):
        df_1["mgf_order"] = order_range
        df_1["mgf_path"] = classify_path
        df_1['mgf_path_index'] = df_1['mgf_path'] + "/" + df_1['mgf_order'].astype(str)

        col_need = ['posterior_error_probability', 'global_qvalue', 'peptidoform',
                    'mz_array', 'intensity_array', 'charge', 'usi', 'pepmass', 'cluster_accession']

        # 存在则映射，不存在则返回NaN（新建列场景下，NaN即表示跳过赋值）
        df_1['cluster_accession'] = df_1['mgf_path_index'].apply(
            lambda x: map_dict.get(x, np.nan)
        )
        valid_df = df_1[pd.notna(df_1['cluster_accession'])]  # 等价于 df_1[~pd.isna(df_1['cluster_accession'])]

        # 步骤2：选取需要的列col_need
        result_df = valid_df.loc[:, col_need]

        return result_df

    @staticmethod
    def read_cluster_tsv(tsv_file_path: str):
        """
        read and rename col
        :param tsv_file_path:
        :return:
        """
        clu_df = pd.read_csv(tsv_file_path, sep='\t', header=None)
        clu_df.columns = ['mgf_path', 'index', 'cluster_accession']
        clu_df.dropna(axis=0, inplace=True)  # 删除空行
        return clu_df

    def inject_cluster_info(self, path_parquet, clu_map_dict, path_sdrf, spectra_num=1000000, batch_size=200000):
        # read parquet file and the cluster result tsv file of Maracluster program.
        print("inject cluster_info: path_parquet:", path_parquet)
        if type(path_parquet) == list:
            path_parquet = path_parquet[0]
        parquet_file = pq.ParquetFile(path_parquet)

        # cluster_res_df = self.read_cluster_tsv(path_cluster_tsv)
        cluster_res_lst = []

        sample_info_dict = SdrfUtil.get_metadata_dict_from_sdrf(path_sdrf)
        # TODO: 这里后面的结果文件应该要修改为增加他的一个分类情况路径在聚类结果文件当中，这个部分应该在NextFlow里面解决，这里只是为了方便调试
        # cluster_res_df.loc[:, "mgf_path"] = cluster_res_df.loc[:, "mgf_path"].apply(lambda x: 'Homo sapiens/Q '
        #                                                                                       'Exactive/charge3/mgf '
        #                                                                                       'files/' + x)
        # # cluster_res_df.set_index(['mgf_path', 'index'], inplace=True)
        # cluster_res_df.index = cluster_res_df.apply(lambda row: f"{row['mgf_path']}/{row['index']}", axis=1)
        basename = ParquetPathHandler(path_parquet).get_item_info()  # 'PXD008467'

        write_count_dict = defaultdict(int)  # Counting dictionary
        file_index_dict = defaultdict(int)  # the file index dictionary
        SPECTRA_NUM = spectra_num  # The spectra capacity of one mgf
        BATCH_SIZE = batch_size

        for parquet_batch in parquet_file.iter_batches(batch_size=BATCH_SIZE,
                                                       columns=UseCol.PARQUET_COL_TO_MSP.value +
                                                               UseCol.PARQUET_COL_TO_FILTER.value):
            row_group = parquet_batch.to_pandas()
            row_group.rename({'exp_mass_to_charge': 'pepmass'}, axis=1, inplace=True)

            # spectrum
            mgf_group_df = row_group.loc[:, self.combine_info_col]
            #mgf_group_df['usi'] = row_group.apply(lambda row: Parquet2Mgf.get_usi(row, basename), axis=1)  # usi
            mgf_group_df['usi'] = row_group['USI']

            mgf_group_df['mgf_file_path'] = row_group.apply(
                lambda row: '/'.join(sample_info_dict.get(Parquet2Mgf.get_filename_from_usi(row)) +
                                     ['charge' + str(row["charge"]), 'mgf files']), axis=1)

            for group, group_df in mgf_group_df.groupby('mgf_file_path'):
                base_mgf_path = group
                mgf_file_path = (f"{group}/{Path(path_parquet).parts[-1].split('.')[0]}_"
                                 f"{file_index_dict[base_mgf_path] + 1}.mgf")      # 'Homo sapiens/Orbitrap Fusion Lumos/charge4/mgf files/PXD004732-consensus_1.mgf/1': 1

                if write_count_dict[group] + group_df.shape[0] <= SPECTRA_NUM:
                    mgf_order_range = range(write_count_dict[group],
                                            write_count_dict[group] + group_df.shape[0],
                                            1)

                    # 把所有的合并完的信息的df加入到列表当中
                    cluster_res_df = self.map_strategy_to_cluster(clu_map_dict, group_df, mgf_file_path,
                                                                  mgf_order_range)
                    cluster_res_lst.append(cluster_res_df)

                    write_count_dict[group] += group_df.shape[0]
                else:
                    remain_num = SPECTRA_NUM - write_count_dict[group]  # 一个mgf文件剩余的容量
                    if remain_num > 0:
                        group_df_remain = group_df.head(remain_num)
                        mgf_order_range = range(write_count_dict[group],
                                                write_count_dict[group] + group_df_remain.shape[0], 1)
                        cluster_res_df = self.map_strategy_to_cluster(clu_map_dict, group_df, mgf_file_path,
                                                                      mgf_order_range)
                        cluster_res_lst.append(cluster_res_df)

                        write_count_dict[group] += group_df_remain.shape[0]

                    file_index_dict[base_mgf_path] += 1
                    write_count_dict[group] = 0

                    mgf_file_path = (f"{base_mgf_path}/{Path(path_parquet).parts[-1].split('.')[0]}_"
                                     f"{file_index_dict[base_mgf_path] + 1}.mgf")

                    group_df_tail = group_df.tail(group_df.shape[0] - remain_num)
                    mgf_order_range = range(write_count_dict[group],
                                            write_count_dict[group] + group_df_tail.shape[0], 1)
                    cluster_res_df = self.map_strategy_to_cluster(clu_map_dict, group_df, mgf_file_path,
                                                                  mgf_order_range)
                    cluster_res_lst.append(cluster_res_df)

                    write_count_dict[group] += group_df_tail.shape[0]

            return self.combine_res_lst(cluster_res_lst)

    def combine_res_lst(self, lst: list) -> pd.DataFrame:
        df_num = len(lst)
        if df_num > 1:
            return pd.DataFrame(np.vstack(lst), columns=lst[0].columns)
        else:
            return lst[0]
