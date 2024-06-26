import logging
import pandas as pd
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
        self.combine_info_col = ['mz_array', 'intensity_array', 'charge',
                                 'pepmass', 'posterior_error_probability', 'global_qvalue']

    @staticmethod
    def map_strategy_to_cluster(df: pd.DataFrame, df_1: pd.DataFrame, classify_path: str, order_range: range):
        df_1["mgf_order"] = order_range
        df_1["mgf_path"] = classify_path
        df_1.set_index(['mgf_path', 'mgf_order'], inplace=True)
        # 给聚类结果文件加上spectrum和pep和global列
        df_1.index = df.index.to_flat_index()

        # 然后，你可以使用一个循环来将"global value"列和"pep"列添加到df1中
        for map_col in ['posterior_error_probability', 'global_qvalue',
                        'mz_array', 'intensity_array', 'charge', 'usi', 'pepmass']:
            df[map_col] = df.index.to_flat_index().map(df_1[map_col])
        return df

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

    def inject_cluster_info(self, path_parquet, path_sdrf,  path_cluster_tsv, spectra_num=1000000, batch_size=200000):
        # read parquet file and the cluster result tsv file of Maracluster program.
        parquet_file = pq.ParquetFile(path_parquet)
        cluster_res_df = self.read_cluster_tsv(path_cluster_tsv)

        # TODO: 这里后面的结果文件应该要修改为增加他的一个分类情况路径在聚类结果文件当中，这个部分应该在NextFlow里面解决，这里只是为了方便调试
        # cluster_res_df.loc[:, "mgf_path"] = cluster_res_df.loc[:, "mgf_path"].apply(lambda x: 'Homo sapiens/Q '
        #                                                                                       'Exactive/charge2/mgf '
        #                                                                                       'files/' + x)
        cluster_res_df.set_index(['mgf_path', 'index'], inplace=True)

        basename = ParquetPathHandler(path_parquet).get_item_info()  # 'PXD008467'

        write_count_dict = defaultdict(int)  # Counting dictionary
        file_index_dict = defaultdict(int)  # the file index dictionary
        SPECTRA_NUM = spectra_num  # The spectra capacity of one mgf
        BATCH_SIZE = batch_size

        for parquet_batch in parquet_file.iter_batches(batch_size=BATCH_SIZE,
                                                       columns=UseCol.PARQUET_COL_TO_MGF.value +
                                                               UseCol.PARQUET_COL_TO_FILTER.value):
            row_group = parquet_batch.to_pandas()
            row_group.rename({'exp_mass_to_charge': 'pepmass'}, axis=1, inplace=True)

            # spectrum
            mgf_group_df = row_group.loc[:, self.combine_info_col]
            mgf_group_df['usi'] = row_group.apply(lambda row: Parquet2Mgf.get_usi(row, basename), axis=1)  # usi
            # =============================如果使用平均质谱方法来生成共识谱，需要保留时间这一列
            # mgf_group_df['retention_time'] = row_group['retention_time']

            sample_info_dict = SdrfUtil.get_metadata_dict_from_sdrf(path_sdrf)

            mgf_group_df['mgf_file_path'] = row_group.apply(
                lambda row: '/'.join(sample_info_dict.get(row['reference_file_name']) +
                                     ['charge' + str(row["charge"]), 'mgf files']), axis=1)

            for group, group_df in mgf_group_df.groupby('mgf_file_path'):
                base_mgf_path = group
                mgf_file_path = (f"{group}/{Path(path_parquet).parts[-1].split('.')[0]}_"
                                 f"{file_index_dict[base_mgf_path] + 1}.mgf")

                if write_count_dict[group] + group_df.shape[0] <= SPECTRA_NUM:
                    mgf_order_range = range(write_count_dict[group],
                                            write_count_dict[group] + group_df.shape[0],
                                            1)
                    cluster_res_df = self.map_strategy_to_cluster(cluster_res_df, group_df, mgf_file_path, mgf_order_range)

                    write_count_dict[group] += group_df.shape[0]
                else:
                    remain_num = SPECTRA_NUM - write_count_dict[group]  # 一个mgf文件剩余的容量
                    if remain_num > 0:
                        group_df_remain = group_df.head(remain_num)
                        mgf_order_range = range(write_count_dict[group],
                                                write_count_dict[group] + group_df_remain.shape[0], 1)
                        cluster_res_df = self.map_strategy_to_cluster(cluster_res_df, group_df, mgf_file_path, mgf_order_range)

                        write_count_dict[group] += group_df_remain.shape[0]

                    file_index_dict[base_mgf_path] += 1
                    write_count_dict[group] = 0

                    mgf_file_path = (f"{base_mgf_path}/{Path(path_parquet).parts[-1].split('.')[0]}_"
                                     f"{file_index_dict[base_mgf_path] + 1}.mgf")

                    group_df_tail = group_df.tail(group_df.shape[0] - remain_num)
                    mgf_order_range = range(write_count_dict[group],
                                            write_count_dict[group] + group_df_tail.shape[0], 1)
                    cluster_res_df = self.map_strategy_to_cluster(cluster_res_df, group_df, mgf_file_path, mgf_order_range)

                    write_count_dict[group] += group_df_tail.shape[0]

            return cluster_res_df


