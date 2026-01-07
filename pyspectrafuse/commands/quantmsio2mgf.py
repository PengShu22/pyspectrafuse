from multiprocessing import Pool
from typing import Tuple
import click
from pyspectrafuse.mgf_convert.parquet2mgf import Parquet2Mgf
from pyspectrafuse.common.parquet_utils import ParquetPathHandler
from pyspectrafuse.common.sdrf_utils import SdrfUtil
import logging
from pathlib import Path
import os

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

logger = logging.getLogger(__name__)


@click.command("convert-mgf", short_help="Convert parquet files to MGF format")
@click.option('--parquet_dir', '-p', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), help='The directory where the parquet files are located')
# @click.option('--sdrf_file_path', '-s', help='The path to the sdrf file')
# @click.option('--output_path', '-o', help='The output directory')
@click.option('--batch_size', '-b', default=100000, help='The batch size of each parquet pass')
@click.option('--spectra_capacity', '-c', default=1000000, help='Number of spectra on each MGF file')
@click.option('--task_parallel', '-t', default=1, help='The number of parquet files that can be converted in parallel')
def quantmsio2mgf(parquet_dir: str, batch_size: int = 100000, 
                  spectra_capacity: int = 1000000, task_parallel: int = 1) -> None:
    """Convert all parquet files in the specified directory to MGF format.
    
    The conversion is based on the SDRF file and the original parquet files
    from the experiment. Files are processed in parallel using multiprocessing.
    
    Args:
        parquet_dir: Directory containing parquet files
        batch_size: Batch size for reading parquet files (default: 100000)
        spectra_capacity: Maximum number of spectra per MGF file (default: 1000000)
        task_parallel: Number of parallel tasks (default: 1)
        
    Raises:
        FileNotFoundError: If no parquet files or SDRF file found
    """
    parquet_file_path_lst = ParquetPathHandler.iter_parquet_dir(parquet_dir)
    if not parquet_file_path_lst:
        raise FileNotFoundError(f"No parquet files found in {parquet_dir}")
    
    sdrf_file_path = SdrfUtil.get_sdrf_file_path(parquet_dir)
    res_file_path = str(Path(parquet_dir) / 'mgf_output')

    pool = Pool(processes=os.cpu_count())  # os.cpu_count()

    # Pack all required parameters into a list of tuples
    tasks = [(parquet_file_path, sdrf_file_path, res_file_path, batch_size, spectra_capacity) for parquet_file_path
             in
             parquet_file_path_lst]

    # Use pool.imap to execute tasks in parallel with specified maximum parallel tasks
    for _ in pool.imap(Parquet2Mgf().convert_to_mgf_task, tasks, chunksize=task_parallel):
        pass

    # Close the process pool and wait for all processes to complete
    pool.close()
    pool.join()

    logger.info("All tasks have completed")
