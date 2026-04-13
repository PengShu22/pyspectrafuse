from typing import Tuple
import pyarrow.parquet as pq
import pandas as pd
import numpy as np
import ast
from pathlib import Path
from collections import defaultdict
import logging
from pyspectrafuse.common.constant import UseCol, ParquetSchemaAdapter
import os
from pyspectrafuse.common.parquet_utils import ParquetPathHandler

logger = logging.getLogger(__name__)


class Parquet2Mgf:

    @staticmethod
    def write2mgf(target_path: str, write_content: str):
        # NOTE: This function may be called multiple times for the same target file (parquet batches).
        # Always ensure we separate batches with blank lines, otherwise we can end up with:
        #   END IONSBEGIN IONS
        # which breaks strict MGF parsers (e.g. MaRaCluster).
        needs_separator = os.path.exists(target_path) and os.path.getsize(target_path) > 0
        with open(target_path, 'a') as f:
            logger.info(f"Writing spectrum to MGF file at path: {target_path}")
            if needs_separator:
                f.write("\n")
            f.write(write_content)
            # Ensure file ends with a newline to make subsequent appends safe.
            if write_content and not write_content.endswith("\n"):
                f.write("\n")

    @staticmethod
    def get_mz_intensity_str(mz_series, intensity_series) -> str:
        """
        Combine the m/z and intensity arrays into a single string
        :param mz_series: m/z array(string type)
        :param intensity_series: intensity array(string type)
        :return: Combined string of m/z and intensity
        """
        if mz_series is None or intensity_series is None:
            return ""

        # Parse both strings to lists at once to avoid repeated parsing
        mz_list = mz_series
        intensity_list = intensity_series

        combined_list = [f"{mz} {intensity}" for mz, intensity in zip(mz_list, intensity_list)]

        # Use the join method to concatenate the list of strings into a single string
        combined_str = '\n'.join(combined_list)
        return combined_str

    @staticmethod
    def get_usi(row: pd.Series, dataset_id: str) -> str:
        """Generate USI string from row data.

        Args:
            row: Pandas Series with spectrum data
            dataset_id: Dataset identifier (e.g., PXD008467)

        Returns:
            USI string
        """
        usi_str = (f'mzspec:{dataset_id}:{row["reference_file_name"]}:'
                   f'scan:{str(row["scan"])}:{row["sequence"]}/{row["precursor_charge"]}')
        return usi_str

    @staticmethod
    def get_spectrum(row: pd.Series, dataset_id: str) -> str:
        res_str = (f"BEGIN IONS\n"  # begin
                   f'TITLE=id={row["USI"]}\n'  # usi
                   f'PEPMASS={str(row["exp_mass_to_charge"])}\n'  # pepmass
                   f'CHARGE={str(row["charge"])}+\n'  # charge
                   f'{Parquet2Mgf.get_mz_intensity_str(row["mz_array"], row["intensity_array"])}\n'  # mz and intensity
                   f'END IONS'  # end
                   )

        return res_str

    @staticmethod
    def get_filename_from_usi(row: pd.Series) -> str:
        """Extract filename from USI string.

        Args:
            row: Pandas Series containing 'USI' column

        Returns:
            Filename extracted from USI
        """
        return row['USI'].split(':')[2]

    def convert_to_mgf(self, parquet_path: str, sample_info_dict: dict, output_path: str,
                       batch_size: int, spectra_capacity: int,
                       skip_instrument: bool = False) -> None:
        """Convert a parquet file to MGF format, grouped by species/instrument/charge.

        Args:
            parquet_path: Path to the PSM parquet file.
            sample_info_dict: {run_file_name: [species, instrument]} metadata dict.
            output_path: Output directory for MGF files.
            batch_size: Batch size for reading parquet (default 100000).
            spectra_capacity: Max spectra per MGF file (default 1000000).
            skip_instrument: If True, all runs use 'all_instruments'.
        """
        Path(output_path).mkdir(parents=True, exist_ok=True)
        # loading Parquet file
        parquet_file = pq.ParquetFile(parquet_path)

        write_count_dict = defaultdict(int)  # Counting dictionary
        relation_dict = defaultdict(int)  # the file index dictionary
        SPECTRA_NUM = spectra_capacity  # The spectra capacity of one mgf
        BATCH_SIZE = batch_size  # The batch size of each parquet pass

        # Determine which columns to read (supports both MSNet and QPX schemas)
        parquet_columns = [f.name for f in parquet_file.schema_arrow]
        # All candidate column names from both formats
        mgf_candidates = UseCol.PARQUET_COL_TO_MGF.value + ['charge', 'observed_mz', 'run_file_name', 'peptidoform']
        read_cols = ParquetSchemaAdapter.available_columns(parquet_columns, mgf_candidates)

        for parquet_batch in parquet_file.iter_batches(batch_size=BATCH_SIZE,
                                                       columns=read_cols):
            mgf_group_df = pd.DataFrame()
            row_group = parquet_batch.to_pandas()

            # Normalize column names (handles both MSNet and QPX schemas)
            row_group = ParquetSchemaAdapter.adapt(row_group)

            # Ensure canonical column names exist for USI and spectrum generation
            # After adapt(): charge, pepmass (observed_mz), reference_file_name (run_file_name)
            if 'precursor_charge' in row_group.columns and 'charge' not in row_group.columns:
                row_group = row_group.rename(columns={'precursor_charge': 'charge'})
            if 'exp_mass_to_charge' in row_group.columns and 'pepmass' not in row_group.columns:
                row_group = row_group.rename(columns={'exp_mass_to_charge': 'pepmass'})

            # Build USI: need reference_file_name, scan, sequence, charge
            dataset_id = ParquetPathHandler(parquet_path).get_item_info()
            row_group['USI'] = ('mzspec:' + dataset_id + ':' +
                                row_group['reference_file_name'].astype(str) + ':scan:' +
                                row_group['scan'].astype(str) + ':' +
                                row_group['sequence'].astype(str) + '/' +
                                row_group['charge'].astype(str))

            # Rename pepmass to exp_mass_to_charge for get_spectrum compatibility
            if 'pepmass' in row_group.columns and 'exp_mass_to_charge' not in row_group.columns:
                row_group['exp_mass_to_charge'] = row_group['pepmass']

            row_group = row_group.loc[:,
                        ['USI', 'sequence', 'mz_array', 'intensity_array', 'charge', 'exp_mass_to_charge']]

            mgf_group_df['spectrum'] = row_group.apply(
                lambda row: Parquet2Mgf.get_spectrum(row, ParquetPathHandler(parquet_path).get_item_info()), axis=1)

            def _get_mgf_path(row):
                fname = Parquet2Mgf.get_filename_from_usi(row)
                info = sample_info_dict.get(fname)
                if info is None:
                    logger.warning(f"No metadata for run file: {fname}")
                    info = ['Unknown', 'Unknown']
                return '/'.join(info + ['charge' + str(row["charge"]), 'mgf files'])

            mgf_group_df['mgf_file_path'] = row_group.apply(_get_mgf_path, axis=1)

            # Pre-compute basename to avoid repeated string operations
            basename_parquet = Path(parquet_path).parts[-1].split('.')[0]

            for group, group_df in mgf_group_df.groupby('mgf_file_path'):
                base_mgf_path = f"{output_path}/{group}"
                file_index = relation_dict[base_mgf_path] + 1
                mgf_file_path = f"{base_mgf_path}/{basename_parquet}_{file_index}.mgf"
                Path(mgf_file_path).parent.mkdir(parents=True, exist_ok=True)

                if write_count_dict[group] + group_df.shape[0] <= SPECTRA_NUM:
                    Parquet2Mgf.write2mgf(mgf_file_path, '\n\n'.join(group_df["spectrum"]))
                    write_count_dict[group] += group_df.shape[0]
                else:
                    remain_num = SPECTRA_NUM - write_count_dict[group]
                    if remain_num > 0:
                        group_df_remain = group_df.head(remain_num)
                        Parquet2Mgf.write2mgf(mgf_file_path, '\n\n'.join(group_df_remain["spectrum"]))
                        write_count_dict[group] += group_df_remain.shape[0]

                    # Update the index of the read and mgf files
                    relation_dict[base_mgf_path] += 1
                    write_count_dict[group] = 0

                    file_index = relation_dict[base_mgf_path] + 1
                    mgf_file_path = f"{base_mgf_path}/{basename_parquet}_{file_index}.mgf"
                    group_df_tail = group_df.tail(group_df.shape[0] - remain_num)
                    Parquet2Mgf.write2mgf(mgf_file_path, '\n\n'.join(group_df_tail["spectrum"]))
                    write_count_dict[group] += group_df_tail.shape[0]

    def convert_to_mgf_task(self, args) -> None:
        """Task wrapper for parallel processing.

        Args:
            args: Tuple of (parquet_file_path, sample_info_dict, res_file_path,
                  batch_size, spectra_capacity) or with optional 6th element
                  skip_instrument (bool).
        """
        if len(args) >= 6:
            parquet_file_path, sample_info_dict, res_file_path, batch_size, spectra_capacity, skip_instrument = args[:6]
        else:
            parquet_file_path, sample_info_dict, res_file_path, batch_size, spectra_capacity = args
            skip_instrument = False
        logger.info(f"Converting {os.path.basename(parquet_file_path)} to MGF format...")
        self.convert_to_mgf(parquet_file_path, sample_info_dict, res_file_path,
                            batch_size, spectra_capacity,
                            skip_instrument=skip_instrument)
