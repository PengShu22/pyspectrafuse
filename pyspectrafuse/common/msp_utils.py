import uuid
import numpy as np
from typing import Union, Dict, Any
import pandas as pd
from pyspectrafuse.mgf_convert.parquet2mgf import Parquet2Mgf
import logging
import ast

logger = logging.getLogger(__name__)


class MspUtil:
    def __init__(self):
        self.namespace: uuid.UUID = uuid.NAMESPACE_URL

    def usi_to_uuid(self, usi_lst: Union[str, list]) -> uuid.UUID:
        """Convert USI string or list to UUID.
        
        Args:
            usi_lst: USI string or list to convert
            
        Returns:
            UUID5 generated from USI
        """
        return uuid.uuid5(self.namespace, usi_lst if isinstance(usi_lst, str) else str(usi_lst))

    @staticmethod
    def get_num_peaks(array: np.ndarray) -> int:
        """Get number of peaks from array.
        
        Args:
            array: NumPy array containing peak data
            
        Returns:
            Number of peaks (array length)
        """
        return array.shape[0]

    @staticmethod
    def get_msp_dict(name: str, mw: float, comment: str, num_peaks: int, 
                     mz_arr: np.ndarray, intensity_arr: np.ndarray) -> Dict[str, Any]:
        """Create MSP format dictionary.
        
        Args:
            name: Spectrum name
            mw: Molecular weight
            comment: Comment string
            num_peaks: Number of peaks
            mz_arr: m/z array
            intensity_arr: Intensity array
            
        Returns:
            Dictionary in MSP format
        """
        msp_dict = {
            'params': {
                'Name': name,
                'MW': mw,
                'Comment': comment,
                'Num peaks': num_peaks
            },
            'mz_array': mz_arr,
            'intensity_array': intensity_arr
        }
        return msp_dict

    @staticmethod
    def get_msp_fmt(row: pd.Series) -> str:
        """Generate MSP format string from row data.
        
        Args:
            row: Pandas Series containing spectrum data
            
        Returns:
            MSP format string
        """
        # if strategy_type == 'most' or strategy_type == 'best':
        #     name_val = row['usi'].split(':')[-1]
        # elif strategy_type == 'average' or strategy_type == 'bin':
        #     name_val = name_val = ';'.join([i.split(':')[-1] for i in row['usi'].split(';')])
        name_val = row['peptidoform']
        mw_val = row['pepmass']
        num_peaks_val = len(ast.literal_eval(row['mz_array']))
        comment_val = f'clusterID={MspUtil().usi_to_uuid(row["usi"])} Nreps={row["Nreps"]} PEP={row["posterior_error_probability"]}'
        mz_intensity_val = Parquet2Mgf.get_mz_intensity_str(row['mz_array'], row['intensity_array'])

        msp_str_fmt = (f"Name: {name_val}\n"
                       f"MW: {mw_val}\n"
                       f"Comment: {comment_val}\n"
                       f"Num peaks: {num_peaks_val}\n"
                       f"{mz_intensity_val}")

        return msp_str_fmt

    @staticmethod
    def write2msp(target_path: str, write_content: str) -> None:
        """Write MSP format content to file.
        
        Args:
            target_path: Path to output file
            write_content: Content to write (MSP format string)
        """
        with open(target_path, mode='a') as f:
            logger.info(f'Writing MSP format spectrum to file: {target_path}')
            f.write(write_content)







