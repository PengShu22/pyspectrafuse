import logging
from typing import List, Dict, Tuple
import click
import uuid
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pyspectrafuse.common.constant import UseCol
from pyspectrafuse.common.parquet_utils import ParquetPathHandler
from pyspectrafuse.common.schemas import CLUSTER_META_SCHEMA, PSM_MEMBERSHIP_SCHEMA
from pyspectrafuse.common.duckdb_ops import compute_stats_and_purity
from pyspectrafuse.cluster_parquet_combine.cluster_res_handler import ClusterResHandler
from pyspectrafuse.cluster_parquet_combine.combine_cluster_and_parquet import CombineCluster2Parquet
from pyspectrafuse.commands.spectrum2msp import create_consensus_strategy, find_target_ext_files

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
logger = logging.getLogger(__name__)


def build_cluster_parquet(
    parquet_dir: str,
    cluster_tsv_file: str,
    species: str,
    instrument: str,
    charge: str,
    method_type: str = 'best',
    output_dir: str = None,
    **strategy_kwargs
) -> Tuple[str, str]:
    """Build cluster_metadata and psm_cluster_membership parquet files.

    Args:
        parquet_dir: Project directory containing parquet and SDRF files
        cluster_tsv_file: Path to MaRaCluster output TSV file
        species: Species name
        instrument: Instrument name
        charge: Charge value
        method_type: Consensus spectrum generation method
        output_dir: Output directory (defaults to parquet_dir/cluster_db)
        **strategy_kwargs: Additional strategy-specific parameters

    Returns:
        Tuple of (cluster_metadata_path, psm_membership_path)
    """
    # Find required files
    from pyspectrafuse.common.qpx_metadata import get_metadata_dict
    sample_info_dict = get_metadata_dict(parquet_dir)

    path_parquet_lst = find_target_ext_files(parquet_dir, '.parquet')
    if not path_parquet_lst:
        raise FileNotFoundError(f"No parquet files found in {parquet_dir}")

    # Setup output directory
    if output_dir is None:
        output_dir = str(Path(parquet_dir) / 'cluster_db')
    out_path = Path(output_dir) / species / instrument / charge
    out_path.mkdir(parents=True, exist_ok=True)

    # Load cluster results and combine with parquet
    cluster_res_dict = ClusterResHandler.get_cluster_dict(
        Path(cluster_tsv_file), species, instrument, charge)
    logger.info(f"Loaded cluster dictionary with {len(cluster_res_dict)} entries")

    combiner = CombineCluster2Parquet()
    df = combiner.inject_cluster_info(
        path_parquet=path_parquet_lst,
        sample_info_dict=sample_info_dict,
        clu_map_dict=cluster_res_dict)

    if df.empty:
        logger.warning("No spectra matched cluster assignments")
        return None, None

    # Extract project accession from parquet path
    basename = ParquetPathHandler(path_parquet_lst[0]).get_item_info()

    # --- Build PSM cluster membership table ---
    psm_df = _build_psm_membership(df, species, instrument, basename)

    # --- Build consensus spectra and cluster metadata ---
    consensus_strategy = create_consensus_strategy(method_type=method_type, **strategy_kwargs)
    consensus_df, single_df = consensus_strategy.consensus_spectrum_aggregation(df)

    cluster_meta_df = _build_cluster_metadata(
        consensus_df, single_df, psm_df, species, instrument, charge, method_type)

    # Write parquet files
    meta_path = str(out_path / 'cluster_metadata.parquet')
    psm_path = str(out_path / 'psm_cluster_membership.parquet')

    _write_cluster_metadata_parquet(cluster_meta_df, meta_path)
    _write_psm_membership_parquet(psm_df, psm_path)

    logger.info(f"Cluster metadata written to {meta_path} ({len(cluster_meta_df)} clusters)")
    logger.info(f"PSM membership written to {psm_path} ({len(psm_df)} PSMs)")

    return meta_path, psm_path


def _build_psm_membership(df: pd.DataFrame, species: str, instrument: str,
                           project_accession: str) -> pd.DataFrame:
    """Build PSM-to-cluster membership dataframe (vectorized)."""
    psm_df = pd.DataFrame({
        'cluster_id': df['cluster_accession'].astype(str).values,
        'usi': df['usi'].astype(str).values,
        'project_accession': project_accession,
        'reference_file_name': df['reference_file_name'].astype(str).values if 'reference_file_name' in df.columns else '',
        'scan': pd.to_numeric(df['scan'], errors='coerce').fillna(0).astype(int).values if 'scan' in df.columns else 0,
        'peptidoform': df['peptidoform'].astype(str).values,
        'charge': pd.to_numeric(df['charge'], errors='coerce').fillna(0).astype(int).values,
        'precursor_mz': pd.to_numeric(df['pepmass'], errors='coerce').astype(float).values,
        'posterior_error_probability': pd.to_numeric(df['posterior_error_probability'], errors='coerce').astype(float).values,
        'global_qvalue': pd.to_numeric(df['global_qvalue'], errors='coerce').astype(float).values,
        'species': species,
        'instrument': instrument,
    })
    return psm_df


def _build_cluster_metadata(consensus_df: pd.DataFrame, single_df: pd.DataFrame,
                             psm_df: pd.DataFrame, species: str, instrument: str,
                             charge: str, method_type: str) -> pd.DataFrame:
    """Build cluster metadata dataframe with consensus spectra.

    Uses DuckDB for stats/purity aggregation, pandas for spectrum array assembly.
    """
    # Cluster stats + purity via DuckDB
    stats_purity = compute_stats_and_purity(psm_df)

    # Build metadata from consensus + single DataFrames (pandas for list columns)
    charge_int = int(charge.replace('charge', '')) if 'charge' in str(charge) else int(charge)
    parts = []

    for source_df in [consensus_df, single_df]:
        if source_df is None or source_df.empty:
            continue

        part = pd.DataFrame()
        part['cluster_id'] = source_df['cluster_accession'].astype(str).values
        part['species'] = species
        part['instrument'] = instrument
        part['charge'] = charge_int
        part['peptidoform'] = source_df['peptidoform'].astype(str).values
        part['peptide_sequence'] = part['peptidoform'].str.split('/').str[0].str.replace(
            r'\[.*?\]', '', regex=True)

        def _to_f32_list(arr):
            if isinstance(arr, np.ndarray):
                return arr.astype(np.float32).tolist()
            if isinstance(arr, str):
                import ast
                return [float(x) for x in ast.literal_eval(arr)]
            if isinstance(arr, (list, tuple)):
                return [float(x) for x in arr]
            if arr is None:
                logger.warning("Null spectrum array encountered — returning empty list")
                return []
            logger.warning(f"Unexpected spectrum array type {type(arr)} — returning empty list")
            return []

        part['consensus_mz_array'] = source_df['mz_array'].apply(_to_f32_list).values
        part['consensus_intensity_array'] = source_df['intensity_array'].apply(_to_f32_list).values
        part['consensus_method'] = method_type
        part['precursor_mz'] = pd.to_numeric(source_df['pepmass'], errors='coerce').astype(float).values
        parts.append(part)

    if not parts:
        return pd.DataFrame()

    meta_df = pd.concat(parts, ignore_index=True)

    # Merge with DuckDB-computed stats and purity
    if not meta_df.empty and not stats_purity.empty:
        meta_df = meta_df.merge(stats_purity, on='cluster_id', how='left')
        meta_df['member_count'] = meta_df['member_count'].fillna(1).astype(int)
        meta_df['project_count'] = meta_df['project_count'].fillna(1).astype(int)
        meta_df['purity'] = meta_df['purity'].fillna(1.0)

    return meta_df


def _write_cluster_metadata_parquet(df: pd.DataFrame, path: str):
    """Write cluster metadata to parquet with canonical schema."""
    if df.empty:
        logger.warning("Empty cluster metadata, skipping write")
        return
    table = pa.Table.from_pandas(df, schema=CLUSTER_META_SCHEMA, preserve_index=False)
    pq.write_table(table, path, compression='zstd')


def _write_psm_membership_parquet(df: pd.DataFrame, path: str):
    """Write PSM cluster membership to parquet with canonical schema."""
    if df.empty:
        logger.warning("Empty PSM membership, skipping write")
        return
    df = df.sort_values('cluster_id')
    table = pa.Table.from_pandas(df, schema=PSM_MEMBERSHIP_SCHEMA, preserve_index=False)
    pq.write_table(table, path, compression='zstd')


@click.command("cluster-parquet", short_help="Generate cluster Parquet database from clustering results")
@click.option('--parquet_dir', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='Project directory containing parquet and SDRF files')
@click.option('--method_type', required=True, type=click.Choice(['best', 'most', 'bin', 'average']),
              help='Consensus spectrum generation method')
@click.option('--cluster_tsv_file', required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='MaRaCluster output TSV file')
@click.option('--species', required=True, help='Species name')
@click.option('--instrument', required=True, help='Instrument name')
@click.option('--charge', required=True, help='Charge value')
@click.option('--output_dir', default=None, type=click.Path(file_okay=False, dir_okay=True),
              help='Output directory (defaults to parquet_dir/cluster_db)')
@click.option('--sim', default='dot', help='Similarity measure (for most method)')
@click.option('--fragment_mz_tolerance', default=0.02, help='Fragment m/z tolerance (for most method)')
@click.option('--min_mz', default=100, help='Minimum m/z (for bin method)')
@click.option('--max_mz', default=2000, help='Maximum m/z (for bin method)')
@click.option('--bin_size', default=0.02, help='Bin size (for bin method)')
@click.option('--peak_quorum', default=0.25, help='Peak quorum (for bin method)')
@click.option('--edge_case_threshold', default=0.5, help='Edge case threshold (for bin method)')
@click.option('--diff_thresh', default=0.01, help='Diff threshold (for average method)')
@click.option('--dyn_range', default=1000, help='Dynamic range (for average method)')
@click.option('--min_fraction', default=0.5, help='Min fraction (for average method)')
@click.option('--pepmass', type=click.Choice(['naive_average', 'neutral_average', 'lower_median']),
              default='lower_median')
@click.option('--msms_avg', type=click.Choice(['naive', 'weighted']), default='weighted')
def cluster2parquet(parquet_dir: str, method_type: str, cluster_tsv_file: str,
                    species: str, instrument: str, charge: str, output_dir: str,
                    sim: str, fragment_mz_tolerance: float,
                    min_mz: float, max_mz: float, bin_size: float,
                    peak_quorum: float, edge_case_threshold: float,
                    diff_thresh: float, dyn_range: int, min_fraction: float,
                    pepmass: str, msms_avg: str) -> None:
    """Generate Parquet database files from clustering results.

    Creates two parquet files:
    - cluster_metadata.parquet: consensus spectra with quality metrics per cluster
    - psm_cluster_membership.parquet: full PSM-to-cluster mapping with provenance
    """
    strategy_kwargs = dict(
        sim=sim, fragment_mz_tolerance=fragment_mz_tolerance,
        min_mz=min_mz, max_mz=max_mz, bin_size=bin_size,
        peak_quorum=peak_quorum, edge_case_threshold=edge_case_threshold,
        diff_thresh=diff_thresh, dyn_range=dyn_range,
        min_fraction=min_fraction, pepmass=pepmass, msms_avg=msms_avg
    )

    meta_path, psm_path = build_cluster_parquet(
        parquet_dir=parquet_dir,
        cluster_tsv_file=cluster_tsv_file,
        species=species,
        instrument=instrument,
        charge=charge,
        method_type=method_type,
        output_dir=output_dir,
        **strategy_kwargs
    )

    if meta_path:
        click.echo(f"Cluster metadata: {meta_path}")
        click.echo(f"PSM membership: {psm_path}")
    else:
        click.echo("No output generated (no matching spectra)")
