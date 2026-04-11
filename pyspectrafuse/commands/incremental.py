"""CLI commands for incremental clustering.

Provides two sub-commands:
- extract-reps: Extract representative spectra from existing cluster DB to MGF.
- merge-clusters: Resolve re-clustering results and merge into existing DB.
"""
import logging
from pathlib import Path

import click

logger = logging.getLogger(__name__)


@click.group("incremental", short_help="Incremental clustering utilities")
def incremental():
    """Incremental clustering: add new projects without re-clustering everything."""


@incremental.command("extract-reps", short_help="Extract representative spectra to MGF")
@click.option('--cluster_metadata', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='Path to existing cluster_metadata.parquet')
@click.option('--output_mgf', required=True,
              type=click.Path(file_okay=True, dir_okay=False),
              help='Output MGF file path for representative spectra')
@click.option('--max_clusters', default=None, type=int,
              help='Limit number of representatives (for testing)')
def extract_reps(cluster_metadata: str, output_mgf: str, max_clusters: int) -> None:
    """Extract one representative spectrum per cluster to MGF for re-clustering."""
    from pyspectrafuse.incremental.representative_mgf import extract_representatives

    count = extract_representatives(
        cluster_metadata_path=cluster_metadata,
        output_mgf_path=output_mgf,
        max_clusters=max_clusters,
    )
    click.echo(f"Extracted {count} representative spectra to {output_mgf}")


@incremental.command("merge-clusters", short_help="Merge incremental clustering results")
@click.option('--cluster_tsv', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='MaRaCluster output TSV from incremental clustering')
@click.option('--rep_mgf', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='Representative MGF used in the incremental clustering run')
@click.option('--existing_metadata', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='Path to existing cluster_metadata.parquet')
@click.option('--existing_membership', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='Path to existing psm_cluster_membership.parquet')
@click.option('--new_parquet_dir', required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='Directory with new project parquet + SDRF files')
@click.option('--species', required=True, help='Species name')
@click.option('--instrument', required=True, help='Instrument name')
@click.option('--charge', required=True, help='Charge value')
@click.option('--method_type', default='best',
              type=click.Choice(['best', 'most', 'bin', 'average']),
              help='Consensus spectrum method')
@click.option('--output_dir', default=None, type=click.Path(file_okay=False, dir_okay=True),
              help='Output directory (defaults to same as existing files)')
def merge_clusters(
    cluster_tsv: str,
    rep_mgf: str,
    existing_metadata: str,
    existing_membership: str,
    new_parquet_dir: str,
    species: str,
    instrument: str,
    charge: str,
    method_type: str,
    output_dir: str,
) -> None:
    """Resolve incremental clustering and merge into existing cluster DB.

    Workflow:
    1. Parse MaRaCluster TSV and resolve new cluster IDs to existing ones
    2. Inject cluster info into new parquet data
    3. Append new PSMs to membership table
    4. Re-run consensus strategy on updated clusters
    5. Rebuild cluster metadata
    """
    from pyspectrafuse.incremental.resolve_clusters import resolve_incremental_clusters
    from pyspectrafuse.incremental.merge_results import (
        append_psm_membership, rebuild_cluster_metadata)
    from pyspectrafuse.cluster_parquet_combine.combine_cluster_and_parquet import CombineCluster2Parquet
    from pyspectrafuse.cluster_parquet_combine.cluster_res_handler import ClusterResHandler
    from pyspectrafuse.commands.spectrum2msp import create_consensus_strategy, find_target_ext_files

    # Step 1: Resolve cluster IDs
    click.echo("Resolving incremental cluster IDs...")
    id_map, new_spectra_df = resolve_incremental_clusters(cluster_tsv, rep_mgf)
    click.echo(f"  Resolved {len(id_map)} cluster IDs, {len(new_spectra_df)} new spectra")

    # Step 2: Load new project data with cluster assignments
    click.echo("Loading new project data...")
    sdrf_files = find_target_ext_files(new_parquet_dir, '.sdrf.tsv')
    if not sdrf_files:
        raise FileNotFoundError(f"No SDRF file found in {new_parquet_dir}")

    parquet_files = find_target_ext_files(new_parquet_dir, '.parquet')
    if not parquet_files:
        raise FileNotFoundError(f"No parquet files found in {new_parquet_dir}")

    # Build cluster map dict using the resolved IDs
    cluster_res_dict = ClusterResHandler.get_cluster_dict(
        Path(cluster_tsv), species, instrument, charge)
    # Remap values through the id_map
    resolved_dict = {k: id_map.get(str(v), str(v)) for k, v in cluster_res_dict.items()}

    combiner = CombineCluster2Parquet()
    new_df = combiner.inject_cluster_info(
        path_parquet=parquet_files,
        clu_map_dict=resolved_dict,
        path_sdrf=sdrf_files[0])

    if new_df.empty:
        click.echo("No new spectra matched cluster assignments")
        return

    # Step 3: Append to PSM membership
    click.echo("Appending new PSMs to membership table...")
    if output_dir is None:
        output_dir = str(Path(existing_membership).parent)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    out_membership = str(Path(output_dir) / 'psm_cluster_membership.parquet')
    # Build a minimal PSM DataFrame from new_df
    from pyspectrafuse.common.parquet_utils import ParquetPathHandler
    import pandas as pd

    basename = ParquetPathHandler(parquet_files[0]).get_item_info()
    psm_new = pd.DataFrame({
        'cluster_id': new_df['cluster_accession'].astype(str),
        'usi': new_df['usi'].astype(str),
        'project_accession': basename,
        'reference_file_name': new_df['reference_file_name'].astype(str) if 'reference_file_name' in new_df.columns else '',
        'scan': pd.to_numeric(new_df['scan'], errors='coerce').fillna(0).astype(int) if 'scan' in new_df.columns else 0,
        'peptidoform': new_df['peptidoform'].astype(str),
        'charge': pd.to_numeric(new_df['charge'], errors='coerce').fillna(0).astype(int),
        'precursor_mz': pd.to_numeric(new_df['pepmass'], errors='coerce').astype(float),
        'posterior_error_probability': pd.to_numeric(new_df['posterior_error_probability'], errors='coerce').astype(float),
        'global_qvalue': pd.to_numeric(new_df['global_qvalue'], errors='coerce').astype(float),
        'species': species,
        'instrument': instrument,
    })

    append_psm_membership(existing_membership, psm_new, out_membership)

    # Step 4: Re-run consensus strategy on all data for affected clusters
    click.echo("Generating consensus spectra...")
    # Read the full merged membership to get all PSMs for affected clusters
    full_membership = pd.read_parquet(out_membership)
    affected_clusters = set(new_df['cluster_accession'].astype(str).unique())

    # For affected clusters, we need to re-run consensus on their full data
    # For unaffected clusters, keep existing metadata
    consensus_strategy = create_consensus_strategy(method_type=method_type)
    consensus_df, single_df = consensus_strategy.consensus_spectrum_aggregation(
        full_membership[full_membership['cluster_id'].isin(affected_clusters)]
        if len(affected_clusters) < len(full_membership['cluster_id'].unique())
        else full_membership
    )

    # Step 5: Rebuild cluster metadata
    click.echo("Rebuilding cluster metadata...")
    out_metadata = str(Path(output_dir) / 'cluster_metadata.parquet')
    rebuild_cluster_metadata(
        psm_membership_path=out_membership,
        consensus_df=consensus_df,
        single_df=single_df,
        species=species,
        instrument=instrument,
        charge=charge,
        method_type=method_type,
        output_path=out_metadata,
    )

    click.echo(f"Incremental merge complete:")
    click.echo(f"  Metadata: {out_metadata}")
    click.echo(f"  Membership: {out_membership}")
