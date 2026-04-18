"""CLI commands for incremental clustering.

Only one sub-command remains: ``extract-reps-dat`` writes one representative
spectrum per cluster from an existing ``cluster_metadata.parquet`` into
MaRaCluster's ``.dat`` format so it can be re-clustered alongside new data.
The merge step is folded into ``build-cluster-db`` — pass it
``--existing_metadata`` and ``--existing_membership`` to merge into an
existing cluster DB instead of building a fresh one.
"""
import logging

import click

logger = logging.getLogger(__name__)


@click.group("incremental", short_help="Incremental clustering utilities")
def incremental():
    """Incremental clustering helpers."""


@incremental.command("extract-reps-dat",
                     short_help="Extract representative spectra to .dat format")
@click.option('--cluster_metadata', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='Path to existing cluster_metadata.parquet')
@click.option('--output_dir', required=True,
              type=click.Path(file_okay=False, dir_okay=True),
              help='Output directory for .dat, .scan_info.dat, .scan_titles.txt')
@click.option('--charge_filter', default=None, type=int,
              help='Only extract representatives with this charge')
@click.option('--file_idx', default=0, type=int,
              help='File index for .dat struct (default 0)')
@click.option('--max_clusters', default=None, type=int,
              help='Limit number of representatives (for testing)')
def extract_reps_dat(cluster_metadata: str, output_dir: str,
                     charge_filter: int, file_idx: int,
                     max_clusters: int) -> None:
    """Extract one representative spectrum per cluster to .dat binary format.

    Writes MaRaCluster-compatible .dat, .scan_info.dat, and .scan_titles.txt
    files. Representative scan titles use the format ``rep:{cluster_id}`` for
    downstream identification during cluster resolution.
    """
    from pyspectrafuse.incremental.representative_dat import extract_representatives_dat

    written, skipped, dat_path, scaninfo_path, titles_path = extract_representatives_dat(
        cluster_metadata_path=cluster_metadata,
        output_dir=output_dir,
        charge_filter=charge_filter,
        file_idx=file_idx,
        max_clusters=max_clusters,
    )
    click.echo(f"Extracted {written} representatives ({skipped} skipped)")
    click.echo(f"  .dat:        {dat_path}")
    click.echo(f"  .scan_info:  {scaninfo_path}")
    click.echo(f"  .scan_titles: {titles_path}")
