from pathlib import Path
from typing import List
import re


class ParquetPathHandler:
    def __init__(self, path: str):
        self.parquet_path: str = path
        self.path_obj: Path = Path(path)

    def get_mgf_filename(self, mgf_file_index: int = 1) -> str:
        """Generate MGF filename from parquet path and index.

        Args:
            mgf_file_index: Index number for the MGF file

        Returns:
            Generated MGF filename string
        """
        filename = f"{self.path_obj.parts[-1].split('.')[0]}_{str(mgf_file_index)}.mgf"
        return filename

    def get_item_info(self) -> str:
        """Extract item identifier (e.g., PXD008467) from parquet path.

        Returns:
            Item identifier (letters + digits prefix from the last path part)
        """
        # 获取路径的最后一个部分（文件名或目录名）
        last_path_part = self.path_obj.parts[-1]

        # 正则匹配：开头的「一个或多个字母 + 一个或多个数字」组合
        # 自动忽略后续的连字符、UUID、扩展名等所有内容
        match = re.search(r'^([A-Za-z]+\d+)', last_path_part)

        if not match:
            raise ValueError(f"Could not extract project ID from path segment: {last_path_part}"
                             f" (Please ensure the path contains a project identifier starting with letters followed by digits.")

        return match.group(1)

    @staticmethod
    def iter_parquet_dir(dir_path: str) -> List[Path]:
        """
        Extract the path information for all parquet files from the parquet Files subdirectory and return a list
        :param dir_path: parquet file's path
        :return:
        """
        parquet_path_lst = []
        directory_path = Path(dir_path)
        parquet_files = directory_path.rglob('*.parquet')

        # Iterate over all matching.parquet files
        for parquet_file in parquet_files:
            parquet_path_lst.append(parquet_file)

        return parquet_path_lst
