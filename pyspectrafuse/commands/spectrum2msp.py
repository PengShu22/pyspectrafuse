import logging
import click
from pyspectrafuse.common.msp_utils import MspUtil
from pyspectrafuse.cluster_parquet_combine.combine_cluster_and_parquet import CombineCluster2Parquet
from pyspectrafuse.consensus_strategy.average_spectrum_strategy import AverageSpectrumStrategy
from pyspectrafuse.consensus_strategy.binning_strategy import BinningStrategy
from pyspectrafuse.consensus_strategy.best_spetrum_strategy import BestSpectrumStrategy
from pyspectrafuse.consensus_strategy.most_similar_strategy import MostSimilarStrategy

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
REVISION = "0.1.1"

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


@click.command("msp", short_help="get msp format file")
@click.option('--path_parquet', help='parquet file path in the project')
@click.option('--path_sdrf', help='sdrf file path for the project')
@click.option('--path_cluster_tsv', help='Cluster result file of maracluster program (.tsv)')
@click.option('--output', help='The output file path in msp format')
@click.option('--strategy_type', default='best', help='Consensus Spectrum generation method')
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
def spectrum2msp(path_parquet, path_sdrf, path_cluster_tsv, strategy_type, output,
                 sim='dot', fragment_mz_tolerance=0.02,  # most method params
                 min_mz=100., max_mz=2000., bin_size=0.02, peak_quorum=0.25, edge_case_threshold=0.5,  # bin
                 diff_thresh=0.01, dyn_range=1000, min_fraction=0.5, pepmass='lower_median', msms_avg='weighted'):
    df = CombineCluster2Parquet().inject_cluster_info(path_parquet=path_parquet,
                                                      path_sdrf=path_sdrf,
                                                      path_cluster_tsv=path_cluster_tsv)

    consensus_strategy = None
    if strategy_type == "best":
        consensus_strategy = BestSpectrumStrategy()

    elif strategy_type == "most":
        consensus_strategy = MostSimilarStrategy(sim=sim, fragment_mz_tolerance=fragment_mz_tolerance)

    elif strategy_type == 'bin':
        consensus_strategy = BinningStrategy(min_mz=min_mz,
                                             max_mz=max_mz,
                                             bin_size=bin_size,
                                             peak_quorum=peak_quorum,
                                             edge_case_threshold=edge_case_threshold)

    elif strategy_type == 'average':
        consensus_strategy = AverageSpectrumStrategy(DIFF_THRESH=diff_thresh,
                                                     DYN_RANGE=dyn_range,
                                                     MIN_FRACTION=min_fraction,
                                                     pepmass=pepmass,
                                                     msms_avg=msms_avg)

    if consensus_strategy:
        consensus_spectrum_df, single_spectrum_df = consensus_strategy.consensus_spectrum_aggregation(df)  # 获得共识谱的df

        # 转化为msp文件需要的格式
        for spectrum_df in [consensus_spectrum_df, single_spectrum_df]:
            spectrum_df['msp_fmt'] = spectrum_df.apply(lambda row:
                                                       MspUtil.get_msp_fmt(row, strategy_type=strategy_type), axis=1)
            MspUtil.write2msp(output, '\n\n'.join(spectrum_df['msp_fmt']))
    else:
        raise ValueError("Unknown strategy type, The current type can only be one of [best, most, bin, average]")


if __name__ == '__main__':
    consensus_strategy_type = 'best'
    msp_output = f'G:/graduation_project/generate-spectrum-library/PXD008467/msp_test/PXD008467_{consensus_strategy_type}_new.msp.txt'
    parquet_file_path = r"G:\graduation_project\generate-spectrum-library\PXD008467\parquet_files\PXD008467-2.parquet"
    sdrf_file_path = 'G:\\graduation_project\\generate-spectrum-library\\PXD008467\\PXD008467.sdrf.tsv'
    res_file_path = 'G:\\graduation_project\\generate-spectrum-library\\PXD008467/mgf_output'
    cluster_res_tsv = (r"G:\graduation_project\generate-spectrum-library\project_cluster\PXD008467\mgf_output\Homo "
                       r"sapiens\Q Exactive\charge2\maracluster_output\MaRaCluster.clusters_p30.tsv")  # 聚类的结果tsv文件

    spectrum2msp(path_parquet=parquet_file_path,
                 path_sdrf=sdrf_file_path,
                 path_cluster_tsv=cluster_res_tsv,
                 strategy_type=consensus_strategy_type,
                 output=msp_output)
