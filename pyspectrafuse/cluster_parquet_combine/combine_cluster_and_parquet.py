from typing import Union, List, Dict
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

logger = logging.getLogger(__name__)


class CombineCluster2Parquet:

    def __init__(self):
        self.combine_info_col = ['mz_array', 'intensity_array', 'charge', 'peptidoform',
                                 'pepmass', 'posterior_error_probability', 'global_qvalue']

    @staticmethod
    def map_strategy_to_cluster(map_dict: Dict[str, str], df_1: pd.DataFrame, 
                                classify_path: str, order_range: range) -> pd.DataFrame:
        df_1["mgf_order"] = order_range
        df_1["mgf_path"] = classify_path
        df_1['mgf_path_index'] = df_1['mgf_path'] + "/" + df_1['mgf_order'].astype(str)

        col_need = ['posterior_error_probability', 'global_qvalue', 'peptidoform',
                    'mz_array', 'intensity_array', 'charge', 'usi', 'pepmass', 'cluster_accession']

        # Map if exists, return NaN if not (NaN means skip assignment for new column)
        # Use vectorized map() instead of apply() for 10-50x speedup
        df_1['cluster_accession'] = df_1['mgf_path_index'].map(map_dict)
        valid_df = df_1[pd.notna(df_1['cluster_accession'])]  # Equivalent to df_1[~pd.isna(df_1['cluster_accession'])]

        # Step 2: Select required columns col_need
        result_df = valid_df.loc[:, col_need]

        return result_df

    @staticmethod
    def read_cluster_tsv(tsv_file_path: str) -> pd.DataFrame:
        """Read and rename columns in cluster TSV file.
        
        Args:
            tsv_file_path: Path to cluster TSV file
            
        Returns:
            DataFrame with renamed columns
        """
        clu_df = pd.read_csv(tsv_file_path, sep='\t', header=None)
        clu_df.columns = ['mgf_path', 'index', 'cluster_accession']
        clu_df.dropna(axis=0, inplace=True)  # Remove empty rows
        return clu_df

    def inject_cluster_info(self, path_parquet: Union[str, List[str]], 
                           clu_map_dict: Dict[str, str], path_sdrf: str,
                           spectra_num: int = 1000000, batch_size: int = 200000) -> pd.DataFrame:
        # read parquet file and the cluster result tsv file of Maracluster program.
        logger.info(f"Injecting cluster info: path_parquet: {path_parquet}")
        if isinstance(path_parquet, list):
            if not path_parquet:
                raise ValueError("path_parquet list is empty")
            path_parquet = path_parquet[0]
        parquet_file = pq.ParquetFile(path_parquet)

        # cluster_res_df = self.read_cluster_tsv(path_cluster_tsv)
        cluster_res_lst = []

        sample_info_dict = SdrfUtil.get_metadata_dict_from_sdrf(path_sdrf)
        basename = ParquetPathHandler(path_parquet).get_item_info()

        write_count_dict = defaultdict(int)  # Counting dictionary
        file_index_dict = defaultdict(int)  # the file index dictionary
        SPECTRA_NUM = spectra_num  # The spectra capacity of one mgf
        BATCH_SIZE = batch_size

        for parquet_batch in parquet_file.iter_batches(batch_size=BATCH_SIZE,
                                                       columns=UseCol.PARQUET_COL_TO_MSP.value +
                                                               UseCol.PARQUET_COL_TO_FILTER.value):
            row_group = parquet_batch.to_pandas()
            row_group.rename({'exp_mass_to_charge': 'pepmass'}, axis=1, inplace=True)

            # Extract spectrum data
            mgf_group_df = row_group.loc[:, self.combine_info_col]
            mgf_group_df['usi'] = row_group['USI']

            # Vectorized string operations - extract filename from USI and map to sample info
            # USI format: mzspec:dataset:filename:scan:sequence/charge
            filenames = row_group['USI'].str.split(':').str[2]
            sample_info = filenames.map(sample_info_dict)
            charges = 'charge' + row_group['charge'].astype(str)
            mgf_group_df['mgf_file_path'] = (sample_info + '/' + charges + '/mgf files')

            # Pre-compute basename to avoid repeated string operations
            basename_parquet = Path(path_parquet).parts[-1].split('.')[0]
            
            for group, group_df in mgf_group_df.groupby('mgf_file_path'):
                base_mgf_path = group
                file_index = file_index_dict[base_mgf_path] + 1
                mgf_file_path = f"{group}/{basename_parquet}_{file_index}.mgf"

                if write_count_dict[group] + group_df.shape[0] <= SPECTRA_NUM:
                    mgf_order_range = range(write_count_dict[group],
                                            write_count_dict[group] + group_df.shape[0],
                                            1)

                    # Add all merged dataframes to the list
                    cluster_res_df = self.map_strategy_to_cluster(clu_map_dict, group_df, mgf_file_path,
                                                                  mgf_order_range)
                    cluster_res_lst.append(cluster_res_df)

                    write_count_dict[group] += group_df.shape[0]
                else:
                    remain_num = SPECTRA_NUM - write_count_dict[group]  # Remaining capacity of one MGF file
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
        if not lst:
            return pd.DataFrame()
        df_num = len(lst)
        if df_num > 1:
            return pd.DataFrame(np.vstack(lst), columns=lst[0].columns)
        else:
            return lst[0]
