from enum import Enum


class UseCol(Enum):
    PARQUET_COL_TO_MGF = ['USI', 'sequence',
                          'mz_array', 'intensity_array', 'charge', 'exp_mass_to_charge']

    PARQUET_COL_TO_MSP = ['USI', 'sequence', 'peptidoform',
                          'mz_array', 'intensity_array', 'charge', 'exp_mass_to_charge']

    PARQUET_COL_TO_FILTER = ['posterior_error_probability', 'global_qvalue']

    PARQUET_COL_TO_USI = ['reference_file_name', 'scan_number', 'sequence', 'charge']

