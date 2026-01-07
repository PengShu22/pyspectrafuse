import logging
from typing import List, Optional
import click
from pathlib import Path
from pyspectrafuse.common.msp_utils import MspUtil
from pyspectrafuse.common.parquet_utils import ParquetPathHandler
from pyspectrafuse.cluster_parquet_combine.combine_cluster_and_parquet import CombineCluster2Parquet
from pyspectrafuse.consensus_strategy.average_spectrum_strategy import AverageSpectrumStrategy
from pyspectrafuse.consensus_strategy.binning_strategy import BinningStrategy
from pyspectrafuse.consensus_strategy.best_spetrum_strategy import BestSpectrumStrategy
from pyspectrafuse.consensus_strategy.most_similar_strategy import MostSimilarStrategy
from pyspectrafuse.cluster_parquet_combine.cluster_res_handler import ClusterResHandler
import uuid

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

logger = logging.getLogger(__name__)


def find_target_ext_files(directory: str, extensions: str) -> List[str]:
    """Find files with specified extension in directory tree.
    
    Args:
        directory: Root directory to search
        extensions: File extension to match (e.g., '.parquet')
        
    Returns:
        List of file paths matching the extension (excluding 'psm' files)
    """
    path = Path(directory)
    res = [str(file) for file in path.rglob(f'*{extensions}') 
           if file.is_file() and 'psm' not in file.name]
    return res


def create_consensus_strategy(method_type: str, sim: str = 'dot', 
                             fragment_mz_tolerance: float = 0.02,
                             min_mz: float = 100., max_mz: float = 2000., 
                             bin_size: float = 0.02, peak_quorum: float = 0.25,
                             edge_case_threshold: float = 0.5,
                             diff_thresh: float = 0.01, dyn_range: int = 1000,
                             min_fraction: float = 0.5, pepmass: str = 'lower_median',
                             msms_avg: str = 'weighted'):
    """Create consensus strategy instance based on method type.
    
    Args:
        method_type: Strategy type ('best', 'most', 'bin', 'average')
        sim: Similarity measure method (for 'most')
        fragment_mz_tolerance: Fragment m/z tolerance (for 'most')
        min_mz: Minimum m/z (for 'bin')
        max_mz: Maximum m/z (for 'bin')
        bin_size: Bin size (for 'bin')
        peak_quorum: Peak quorum threshold (for 'bin')
        edge_case_threshold: Edge case threshold (for 'bin')
        diff_thresh: Difference threshold (for 'average')
        dyn_range: Dynamic range (for 'average')
        min_fraction: Minimum fraction (for 'average')
        pepmass: PEP mass method (for 'average')
        msms_avg: MS/MS averaging method (for 'average')
        
    Returns:
        ConsensusStrategy instance
        
    Raises:
        ValueError: If method_type is not recognized
    """
    if method_type == "best":
        return BestSpectrumStrategy()
    elif method_type == "most":
        return MostSimilarStrategy(sim=sim, fragment_mz_tolerance=fragment_mz_tolerance)
    elif method_type == 'bin':
        return BinningStrategy(min_mz=min_mz, max_mz=max_mz, bin_size=bin_size,
                              peak_quorum=peak_quorum, edge_case_threshold=edge_case_threshold)
    elif method_type == 'average':
        return AverageSpectrumStrategy(DIFF_THRESH=diff_thresh, DYN_RANGE=dyn_range,
                                      MIN_FRACTION=min_fraction, pepmass=pepmass,
                                      msms_avg=msms_avg)
    else:
        raise ValueError(f"Unknown strategy type: {method_type}. "
                         f"Must be one of [best, most, bin, average]")


def write_spectra_to_msp(consensus_df, single_df, output_path: str) -> None:
    """Write consensus and single spectra to MSP file.
    
    Args:
        consensus_df: DataFrame with consensus spectra
        single_df: DataFrame with single spectra
        output_path: Path to output MSP file
    """
    for spectrum_df in [consensus_df, single_df]:
        if spectrum_df.empty:
            logger.warning("Empty spectrum dataframe, skipping")
            continue
        logger.info(f"Processing {len(spectrum_df)} spectra")
        spectrum_df.loc[:, 'msp_fmt'] = spectrum_df.apply(
            lambda row: MspUtil.get_msp_fmt(row), axis=1)
        MspUtil.write2msp(output_path, '\n\n'.join(spectrum_df['msp_fmt']))


@click.command("msp", short_help="get msp format file")
@click.option('--parquet_dir', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), help='The project directory, the directory must obtain parquet and sdrf files.')
@click.option('--method_type', required=True, type=click.Choice(['best', 'most', 'bin', 'average']), help='Consensus Spectrum generation method')
@click.option('--cluster_tsv_file', required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False), help='the MaRaCluster output file')
@click.option('--species', required=True, help='species name')
@click.option('--instrument', required=True, help='instrument name')
@click.option('--charge', required=True, help='charge name')
@click.option('--sim', default='dot', help='The similarity measure method for the most consensus spectrum generation '
                                           'method')
@click.option('--fragment_mz_tolerance', default=0.02,
              help='Fragment m/z tolerance used during spectrum comparison [optional;required for the most_similar method]')
@click.option('--min_mz', default=100,
              help='Minimum m/z to consider for spectrum binning (optional; required for the "bin" method)')
@click.option('--max_mz', default=2000,
              help=' Maximum m/z to consider for spectrum binning (optional; required for the "bin" method).')
@click.option('--bin_size', default=0.02,
              help='Bin size in m/z used for spectrum binning (optional; required for the "bin" method)')
@click.option('--peak_quorum', default=0.25,
              help="Relative number of spectra in a cluster that need to contain a peak for it to be included in the representative spectrum (optional; required for the bin method)")
@click.option('--edge_case_threshold', default=0.5,
              help=' During binning try to correct m/z edge cases where the m/z is closer to the bin edge than the given relative bin size threshold (optional;required for the "bin" method).')
@click.option('--diff_thresh', default=0.01, help='Minimum distance between MS/MS peak clusters.')
@click.option('--dyn_range', default=1000, help='Dynamic range to apply to output spectra')
@click.option('--min_fraction', default=0.5, help='Minimum fraction of cluster spectra where MS/MS peak is present.')
@click.option('--pepmass', type=click.Choice(['naive_average', 'neutral_average', 'lower_median']),
              default='lower_median')
@click.option('--msms_avg', type=click.Choice(['naive', 'weighted']), default='weighted')
def spectrum2msp(parquet_dir: str, method_type: str, cluster_tsv_file: str, 
                 species: str, instrument: str, charge: str,
                 sim: str = 'dot', fragment_mz_tolerance: float = 0.02,
                 min_mz: float = 100., max_mz: float = 2000., bin_size: float = 0.02,
                 peak_quorum: float = 0.25, edge_case_threshold: float = 0.5,
                 diff_thresh: float = 0.01, dyn_range: int = 1000,
                 min_fraction: float = 0.5, pepmass: str = 'lower_median',
                 msms_avg: str = 'weighted') -> None:
    """Generate MSP format file from parquet and cluster data.
    
    This command processes parquet files and cluster results to generate MSP format
    files using various consensus spectrum generation strategies.
    
    Args:
        parquet_dir: Project directory containing parquet and SDRF files
        method_type: Consensus spectrum generation method ('best', 'most', 'bin', 'average')
        cluster_tsv_file: Path to MaRaCluster output TSV file
        species: Species name
        instrument: Instrument name
        charge: Charge value
        sim: Similarity measure method (for 'most' method)
        fragment_mz_tolerance: Fragment m/z tolerance (for 'most' method)
        min_mz: Minimum m/z value (for 'bin' method)
        max_mz: Maximum m/z value (for 'bin' method)
        bin_size: Bin size in m/z (for 'bin' method)
        peak_quorum: Peak quorum threshold (for 'bin' method)
        edge_case_threshold: Edge case threshold (for 'bin' method)
        diff_thresh: Minimum distance between MS/MS peak clusters (for 'average' method)
        dyn_range: Dynamic range to apply (for 'average' method)
        min_fraction: Minimum fraction of cluster spectra (for 'average' method)
        pepmass: PEP mass calculation method (for 'average' method)
        msms_avg: MS/MS averaging method (for 'average' method)
        
    Raises:
        FileNotFoundError: If required files are not found
        ValueError: If method_type is invalid
    """
    # Find required files
    sdrf_files = find_target_ext_files(parquet_dir, '.sdrf.tsv')
    if not sdrf_files:
        raise FileNotFoundError(f"No SDRF file found in {parquet_dir}")
    path_sdrf = sdrf_files[0]
    logger.info(f"SDRF file path: {path_sdrf}")
    
    path_parquet_lst = find_target_ext_files(parquet_dir, '.parquet')
    if not path_parquet_lst:
        raise FileNotFoundError(f"No parquet files found in {parquet_dir}")
    logger.info(f"Found {len(path_parquet_lst)} parquet file(s)")

    # Setup output directory and file
    output_dir = Path(parquet_dir) / 'msp' / species / instrument / charge
    output_dir.mkdir(parents=True, exist_ok=True)
    basename = ParquetPathHandler(parquet_dir).get_item_info()
    output = output_dir / f"{basename}_{uuid.uuid4()}.msp.txt"

    # Load cluster results
    cluster_res_dict = ClusterResHandler.get_cluster_dict(
        Path(cluster_tsv_file), species, instrument, charge)
    logger.info(f"Loaded cluster dictionary with {len(cluster_res_dict)} entries")

    # Combine cluster info with parquet data
    df = CombineCluster2Parquet().inject_cluster_info(
        path_parquet=path_parquet_lst,
        path_sdrf=path_sdrf,
        clu_map_dict=cluster_res_dict)

    # Create consensus strategy and generate spectra
    consensus_strategy = create_consensus_strategy(
        method_type=method_type, sim=sim, fragment_mz_tolerance=fragment_mz_tolerance,
        min_mz=min_mz, max_mz=max_mz, bin_size=bin_size, peak_quorum=peak_quorum,
        edge_case_threshold=edge_case_threshold, diff_thresh=diff_thresh,
        dyn_range=dyn_range, min_fraction=min_fraction, pepmass=pepmass,
        msms_avg=msms_avg)
    
    consensus_spectrum_df, single_spectrum_df = consensus_strategy.consensus_spectrum_aggregation(df)

    # Write spectra to MSP file
    write_spectra_to_msp(consensus_spectrum_df, single_spectrum_df, str(output))
    
    logger.info(f"MSP file generated successfully: {output}")



