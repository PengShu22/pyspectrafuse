from pathlib import Path
from typing import List


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
        """Extract item identifier from parquet path.
        
        Returns:
            Item identifier (typically dataset ID like PXD008467)
        """
        return self.path_obj.parts[-1].split('-')[0]

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
            if parquet_file.parts[-2] == 'parquet_files':
                parquet_path_lst.append(parquet_file)

        return parquet_path_lst
