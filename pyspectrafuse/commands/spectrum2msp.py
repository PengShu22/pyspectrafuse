import logging
from os.path import split

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
REVISION = "0.1.1"

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def find_target_ext_files(directory: str, extensions: str) -> list:
    path = Path(directory)
    res = [str(file) for file in path.rglob(f'*{extensions}') if file.is_file() and 'psm' not in file.name]

    return res


@click.command("msp", short_help="get msp format file")
@click.option('--parquet_dir', help='The project directory, the directory must obtain parquet and sdrf files.')
@click.option('--method_type', help='Consensus Spectrum generation method')
@click.option('--cluster_tsv_file', help='the MaRaCluster output file')
@click.option('--species', help='species name')
@click.option('--instrument', help='instrument name')
@click.option('--charge', help='charge name')
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
def spectrum2msp(parquet_dir, method_type, cluster_tsv_file, species, instrument, charge,
                 sim='dot', fragment_mz_tolerance=0.02,  # most method params
                 min_mz=100., max_mz=2000., bin_size=0.02, peak_quorum=0.25, edge_case_threshold=0.5,  # bin
                 diff_thresh=0.01, dyn_range=1000, min_fraction=0.5, pepmass='lower_median', msms_avg='weighted'):
    """
    :param cluster_tsv_file:
    :param charge:
    :param instrument:
    :param species:
    :param parquet_dir: Project directory name
    :param method_type: Consensus spectrum generation method type
    :param sim:
    :param fragment_mz_tolerance:
    :param min_mz:
    :param max_mz:
    :param bin_size:
    :param peak_quorum:
    :param edge_case_threshold:
    :param diff_thresh:
    :param dyn_range:
    :param min_fraction:
    :param pepmass:
    :param msms_avg:
    :return:
    """
    # print(f'tsv_file_path:{cluster_tsv_file}')
    # print(species)
    # print(instrument)
    # print(charge)
    split_type = ""
    path_sdrf = find_target_ext_files(parquet_dir, '.sdrf.tsv')[0]  # SDRF file path
    print("path_sdrf", path_sdrf)
    path_parquet_lst = find_target_ext_files(parquet_dir, '.parquet')  # parquet_dir is the project path, returns all parquet files
    print("parquet path is",path_parquet_lst)

    output_dir = Path(f'{parquet_dir}/msp/{species}/{instrument}/{charge}')  # Create output MSP directory
    output_dir.mkdir(parents=True, exist_ok=True)
    basename = ParquetPathHandler(parquet_dir).get_item_info()  # PXD008467
    output = f"{output_dir}/{basename}_{uuid.uuid4()}.msp.txt"  # One project corresponds to one MSP format file

    # Get cluster result file, adding species/instrument/charge information since this info is not available in the work directory in Nextflow
    cluster_res_dict = ClusterResHandler.get_cluster_dict(cluster_tsv_file, species, instrument, charge)
    print("cluster_res_dict", cluster_res_dict)

    # TODO: Automatically find parquet files based on charge
    # If parquet files are split by charge, query corresponding parquet files based on charge,
    # as large files may have multiple parquet files per charge
    if split_type == "charge":
        path_parquet = [i for i in path_parquet_lst if i.split('-')[-1][0] == charge[-1]][0]
    else:
        path_parquet = path_parquet_lst

    df = CombineCluster2Parquet().inject_cluster_info(path_parquet=path_parquet, # C:\\E\\graduation\\cluster\\PXD004732\\parquet_files\\PXD004732.parquet
                                                      path_sdrf=path_sdrf, # C:\E\graduation\cluster\PXD004732\PXD004732.sdrf.tsv
                                                      clu_map_dict=cluster_res_dict)  # {'Homo sapiens/Orbitrap Fusion Lumos/charge4/mgf files/PXD004732-consensus_1.mgf/1': 1}

    # Different peptide modifications (based on peptidoform)
    # pep_lst = df['peptidoform'].to_list()
    # print(f"sequence num is {len(pep_lst)}")
    # print(f"unique sequence num is {len(np.unique(pep_lst))}")

    consensus_strategy = None
    if method_type == "best":
        consensus_strategy = BestSpectrumStrategy()

    elif method_type == "most":
        consensus_strategy = MostSimilarStrategy(sim=sim, fragment_mz_tolerance=fragment_mz_tolerance)

    elif method_type == 'bin':
        consensus_strategy = BinningStrategy(min_mz=min_mz,
                                             max_mz=max_mz,
                                             bin_size=bin_size,
                                             peak_quorum=peak_quorum,
                                             edge_case_threshold=edge_case_threshold)

    elif method_type == 'average':
        consensus_strategy = AverageSpectrumStrategy(DIFF_THRESH=diff_thresh,
                                                     DYN_RANGE=dyn_range,
                                                     MIN_FRACTION=min_fraction,
                                                     pepmass=pepmass,
                                                     msms_avg=msms_avg)

    if consensus_strategy:
        consensus_spectrum_df, single_spectrum_df = consensus_strategy.consensus_spectrum_aggregation(
            df)  # Get consensus spectrum dataframe

        # Convert to MSP file format
        for spectrum_df in [consensus_spectrum_df, single_spectrum_df]:
            print("consensus_spectrum_df: ", consensus_spectrum_df)
            print("single_spectrum_df: ", single_spectrum_df)

            spectrum_df.loc[:, 'msp_fmt'] = spectrum_df.apply(lambda row: MspUtil.get_msp_fmt(row), axis=1)
            MspUtil.write2msp(output, '\n\n'.join(spectrum_df['msp_fmt']))
    else:
        raise ValueError("Unknown strategy type, The current type can only be one of [best, most, bin, average]")



