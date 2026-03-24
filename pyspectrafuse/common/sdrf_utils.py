from pathlib import Path
from typing import Dict, List
import pandas as pd
import re
import logging
import numpy as np

logger = logging.getLogger(__name__)


class SdrfUtil:
    def __init__(self, folder_path: str):
        self.folder: str = folder_path

    @staticmethod
    def get_sdrf_file_path(folder: str) -> str:
        """
        get sdrf file path in project folder
        :param folder: folder obtain project files to cluster
        """
        directory_path = Path(folder)
        files = directory_path.rglob('*.sdrf.tsv')
        try:
            return str(next(files))
        except StopIteration:
            raise FileNotFoundError(f'There is no sdrf file in {folder}')

    @staticmethod
    def get_metadata_dict_from_sdrf(sdrf_folder: str, fill_unknown=True):
        """
        从SDRF DataFrame中构建样本信息字典：{数据文件名: [物种, 仪器名称]}

        参数:
            sdrf_feature_df: 包含SDRF标准字段的pandas DataFrame
            fill_unknown: 若为True，空值/无匹配值填充为'Unknown'；若为False，保留np.nan

        返回:
            sample_info_dict: 键为comment[data file]，值为[organism, instrument]的字典
        """

        def extract_instrument(x):
            """内部函数：从comment[instrument]字段提取仪器名称，兼容脏数据"""

            if pd.isna(x) or not isinstance(x, str):
                return np.nan

            match = re.search(r'NT=\s*(.*?)\s*(;|$)', x.strip())
            if not match:
                return np.nan

            return match.group(1).strip()
        sdrf_feature_df = pd.read_csv(sdrf_folder, sep='\t')
        # ---------------------- 1. 清洗并提取仪器名称 ----------------------
        sdrf_feature_df['comment[instrument]'] = sdrf_feature_df['comment[instrument]'].apply(extract_instrument)

        sdrf_feature_df['characteristics[organism]'] = sdrf_feature_df['characteristics[organism]'].apply(
            lambda x: x.strip() if pd.notna(x) and isinstance(x, str) else np.nan
        )

        duplicate_files = sdrf_feature_df[sdrf_feature_df.duplicated('comment[data file]', keep=False)][
            'comment[data file]'].unique()
        if len(duplicate_files) > 0:
            print(
                f"comment[data file] have{len(duplicate_files)} repeat value：{duplicate_files[:3]}")
            # sdrf_feature_df = sdrf_feature_df.drop_duplicates('comment[data file]', keep='first')

        if fill_unknown:
            sdrf_feature_df[['characteristics[organism]', 'comment[instrument]']] = \
                sdrf_feature_df[['characteristics[organism]', 'comment[instrument]']].fillna('Unknown')

        sdrf_feature_df['organism_instrument'] = sdrf_feature_df[
            ['characteristics[organism]', 'comment[instrument]']
        ].values.tolist()
        sdrf_feature_df['comment[data file]'] = sdrf_feature_df['comment[data file]'].apply(lambda x: x.replace('.raw', ''))

        sample_info_dict = sdrf_feature_df.set_index('comment[data file]')['organism_instrument'].to_dict()

        return sample_info_dict


