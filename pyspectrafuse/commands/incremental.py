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
    from pyspectrafuse.incremental.merge_results import append_psm_membership
    from pyspectrafuse.cluster_parquet_combine.combine_cluster_and_parquet import CombineCluster2Parquet
    from pyspectrafuse.cluster_parquet_combine.cluster_res_handler import ClusterResHandler
    from pyspectrafuse.commands.spectrum2msp import find_target_ext_files

    # Step 1: Resolve cluster IDs
    click.echo("Resolving incremental cluster IDs...")
    id_map, new_spectra_df = resolve_incremental_clusters(cluster_tsv, rep_mgf)
    click.echo(f"  Resolved {len(id_map)} cluster IDs, {len(new_spectra_df)} new spectra")

    # Step 2: Load new project data with cluster assignments
    click.echo("Loading new project data...")
    from pyspectrafuse.common.qpx_metadata import get_metadata_dict
    sample_info_dict = get_metadata_dict(new_parquet_dir)

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
        sample_info_dict=sample_info_dict)

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

    # Step 4: Update cluster metadata with new spectra
    click.echo("Updating cluster metadata...")
    import re
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq_write
    from pyspectrafuse.commands.spectrum2msp import create_consensus_strategy

    existing_meta = pd.read_parquet(existing_metadata)
    out_metadata = str(Path(output_dir) / 'cluster_metadata.parquet')

    # Recompute stats from merged membership
    full_membership = pd.read_parquet(out_membership)
    cluster_stats = full_membership.groupby('cluster_id').agg(
        member_count=('usi', 'count'),
        project_count=('project_accession', 'nunique'),
        best_pep=('posterior_error_probability', 'min'),
        best_qvalue=('global_qvalue', 'min'),
    ).reset_index()

    # Purity — double-groupby is 84x faster than lambda x: x.value_counts().iloc[0]
    pair_counts = full_membership.groupby(['cluster_id', 'peptidoform']).size().reset_index(name='_cnt')
    mode_count = pair_counts.groupby('cluster_id')['_cnt'].max()
    total = full_membership.groupby('cluster_id').size()
    purity_series = (mode_count / total)

    # Prepare new spectra grouped by cluster
    new_df_copy = new_df.copy()
    new_df_copy['cluster_accession'] = new_df_copy['cluster_accession'].astype(str)

    existing_meta = existing_meta.copy()
    existing_meta['cluster_id'] = existing_meta['cluster_id'].astype(str)
    existing_ids = set(existing_meta['cluster_id'])

    # Split new spectra: those joining existing clusters vs brand new clusters
    new_for_existing_mask = new_df_copy['cluster_accession'].isin(existing_ids)
    new_for_existing_df = new_df_copy[new_for_existing_mask]
    new_brand_new_df = new_df_copy[~new_for_existing_mask]

    # ── Update existing clusters ──
    existing_clusters_updated = 0
    if len(new_for_existing_df) > 0:
        meta_indexed = existing_meta.set_index('cluster_id')

        if method_type == 'bin':
            # BIN strategy: recompute consensus from existing consensus + new spectra
            # Treat the stored consensus as one additional "spectrum" in the bin
            strategy = create_consensus_strategy('bin')
            affected_ids = new_for_existing_df['cluster_accession'].unique()

            for cid in affected_ids:
                new_spectra = new_for_existing_df[new_for_existing_df['cluster_accession'] == cid]
                if cid not in meta_indexed.index:
                    continue

                existing_row = meta_indexed.loc[cid]
                existing_mz = existing_row['consensus_mz_array']
                existing_int = existing_row['consensus_intensity_array']
                if existing_mz is None or len(existing_mz) == 0:
                    continue

                # Build a mini-DataFrame with existing consensus + new spectra
                existing_mz = np.array(existing_mz, dtype=np.float32)
                existing_int = np.array(existing_int, dtype=np.float32)

                rows = [{
                    'usi': f'existing_consensus:{cid}',
                    'pepmass': float(existing_row['precursor_mz']),
                    'charge': int(existing_row['charge']),
                    'mz_array': existing_mz,
                    'intensity_array': existing_int,
                    'peptidoform': str(existing_row['peptidoform']),
                    'posterior_error_probability': float(existing_row['best_pep']),
                    'global_qvalue': float(existing_row.get('best_qvalue', 1.0)),
                    'cluster_accession': cid,
                }]
                for _, spec in new_spectra.iterrows():
                    pep_str = str(spec['peptidoform'])
                    if '/' not in pep_str:
                        pep_str = pep_str + '/' + str(int(spec['charge']))
                    rows.append({
                        'usi': str(spec['usi']),
                        'pepmass': float(spec['pepmass']),
                        'charge': int(spec['charge']),
                        'mz_array': np.array(spec['mz_array'], dtype=np.float32),
                        'intensity_array': np.array(spec['intensity_array'], dtype=np.float32),
                        'peptidoform': pep_str,
                        'posterior_error_probability': float(spec['posterior_error_probability']),
                        'global_qvalue': float(spec['global_qvalue']),
                        'cluster_accession': cid,
                    })

                mini_df = pd.DataFrame(rows)
                try:
                    consensus_df, _ = strategy.consensus_spectrum_aggregation(mini_df)
                    if not consensus_df.empty:
                        row_out = consensus_df.iloc[0]
                        meta_indexed.loc[cid, 'consensus_mz_array'] = row_out['mz_array'].astype(np.float32).tolist() if isinstance(row_out['mz_array'], np.ndarray) else list(row_out['mz_array'])
                        meta_indexed.loc[cid, 'consensus_intensity_array'] = row_out['intensity_array'].astype(np.float32).tolist() if isinstance(row_out['intensity_array'], np.ndarray) else list(row_out['intensity_array'])
                        meta_indexed.loc[cid, 'precursor_mz'] = float(row_out['pepmass'])
                        existing_clusters_updated += 1
                except Exception as e:
                    logger.warning(f"BIN consensus failed for cluster {cid}: {e}")

            existing_meta = meta_indexed.reset_index().rename(columns={'index': 'cluster_id'})
            if 'cluster_id' not in existing_meta.columns:
                existing_meta = existing_meta.reset_index()

        else:
            # BEST strategy: swap in better-qvalue spectrum if found
            best_new_idx = new_for_existing_df.groupby('cluster_accession')['global_qvalue'].idxmin()
            best_new = new_for_existing_df.loc[best_new_idx].copy()
            best_new = best_new.set_index('cluster_accession')

            meta_indexed = existing_meta.set_index('cluster_id')
            joined = meta_indexed.loc[meta_indexed.index.isin(best_new.index)].copy()
            joined['new_qvalue'] = best_new['global_qvalue'].astype(float)
            better_mask = joined['new_qvalue'] < joined['best_qvalue']
            better_ids = joined[better_mask].index

            if len(better_ids) > 0:
                for cid in better_ids:
                    new_row = best_new.loc[cid]
                    meta_indexed.loc[cid, 'consensus_mz_array'] = new_row['mz_array'].astype(np.float32).tolist() if isinstance(new_row['mz_array'], np.ndarray) else list(new_row['mz_array'])
                    meta_indexed.loc[cid, 'consensus_intensity_array'] = new_row['intensity_array'].astype(np.float32).tolist() if isinstance(new_row['intensity_array'], np.ndarray) else list(new_row['intensity_array'])
                    meta_indexed.loc[cid, 'precursor_mz'] = float(new_row['pepmass'])
                    pep_str = str(new_row['peptidoform'])
                    if '/' not in pep_str:
                        pep_str = pep_str + '/' + str(int(new_row['charge']))
                    meta_indexed.loc[cid, 'peptidoform'] = pep_str
                    meta_indexed.loc[cid, 'peptide_sequence'] = re.sub(r'\[.*?\]', '', pep_str.split('/')[0])
                existing_clusters_updated = len(better_ids)

            existing_meta = meta_indexed.reset_index().rename(columns={'index': 'cluster_id'})
            if 'cluster_id' not in existing_meta.columns:
                existing_meta = existing_meta.reset_index()

    # ── Build metadata for brand-new clusters ──
    new_clusters_added = 0
    if len(new_brand_new_df) > 0:
        charge_int = int(str(charge).replace('charge', ''))
        new_cluster_ids = new_brand_new_df['cluster_accession'].unique()
        new_clusters_added = len(new_cluster_ids)

        if method_type == 'bin' and len(new_brand_new_df) > len(new_cluster_ids):
            # Some new clusters have multiple members — run BIN consensus on them
            strategy = create_consensus_strategy('bin')
            multi_mask = new_brand_new_df['cluster_accession'].isin(
                new_brand_new_df['cluster_accession'].value_counts()[
                    new_brand_new_df['cluster_accession'].value_counts() > 1].index)
            multi_df = new_brand_new_df[multi_mask]
            single_new_df = new_brand_new_df[~multi_mask]

            new_meta_parts = []

            # Process multi-member new clusters with BIN
            if not multi_df.empty:
                # Prepare for consensus strategy
                multi_df = multi_df.copy()
                pep_s = multi_df['peptidoform'].astype(str)
                needs_ch = ~pep_s.str.contains('/', regex=False)
                multi_df['peptidoform'] = pep_s.where(
                    ~needs_ch, pep_s + '/' + multi_df['charge'].astype(int).astype(str))

                try:
                    consensus_df, singles_df = strategy.consensus_spectrum_aggregation(multi_df)
                    for source_df in [consensus_df, singles_df]:
                        if source_df is None or source_df.empty:
                            continue
                        part = pd.DataFrame({
                            'cluster_id': source_df['cluster_accession'].astype(str).values if 'cluster_accession' in source_df.columns else source_df.index.astype(str).values,
                            'species': species,
                            'instrument': instrument,
                            'charge': charge_int,
                            'peptidoform': source_df['peptidoform'].astype(str).values,
                            'peptide_sequence': source_df['peptidoform'].astype(str).str.split('/').str[0].str.replace(r'\[.*?\]', '', regex=True).values,
                            'consensus_mz_array': source_df['mz_array'].values,
                            'consensus_intensity_array': source_df['intensity_array'].values,
                            'consensus_method': method_type,
                            'precursor_mz': pd.to_numeric(source_df['pepmass'], errors='coerce').astype(float).values,
                            'member_count': 1,
                            'project_count': 1,
                            'best_pep': pd.to_numeric(source_df['posterior_error_probability'], errors='coerce').astype(float).values,
                            'best_qvalue': pd.to_numeric(source_df.get('global_qvalue', source_df['posterior_error_probability']), errors='coerce').astype(float).values,
                            'purity': 1.0,
                        })
                        new_meta_parts.append(part)
                except Exception as e:
                    logger.warning(f"BIN consensus failed for new clusters: {e}")
                    # Fall back to best-spectrum for failed clusters
                    single_new_df = pd.concat([single_new_df, multi_df])

            # Process single-member new clusters (same for both strategies)
            if not single_new_df.empty:
                best_idx = single_new_df.groupby('cluster_accession')['global_qvalue'].idxmin()
                best_single = single_new_df.loc[best_idx].set_index('cluster_accession')
                pep_series = best_single['peptidoform'].astype(str)
                needs_charge = ~pep_series.str.contains('/', regex=False)
                pep_series = pep_series.where(
                    ~needs_charge,
                    pep_series + '/' + best_single['charge'].astype(int).astype(str))

                part = pd.DataFrame({
                    'cluster_id': best_single.index.astype(str),
                    'species': species,
                    'instrument': instrument,
                    'charge': charge_int,
                    'peptidoform': pep_series.values,
                    'peptide_sequence': pep_series.str.split('/').str[0].str.replace(
                        r'\[.*?\]', '', regex=True).values,
                    'consensus_mz_array': best_single['mz_array'].values,
                    'consensus_intensity_array': best_single['intensity_array'].values,
                    'consensus_method': method_type,
                    'precursor_mz': best_single['pepmass'].astype(float).values,
                    'member_count': 1,
                    'project_count': 1,
                    'best_pep': best_single['posterior_error_probability'].astype(float).values,
                    'best_qvalue': best_single['global_qvalue'].astype(float).values,
                    'purity': 1.0,
                })
                new_meta_parts.append(part)

            if new_meta_parts:
                new_meta_df = pd.concat(new_meta_parts, ignore_index=True)
                meta = pd.concat([existing_meta, new_meta_df], ignore_index=True)
            else:
                meta = existing_meta
        else:
            # BEST strategy or all single-member: select best spectrum per cluster
            best_new_idx = new_brand_new_df.groupby('cluster_accession')['global_qvalue'].idxmin()
            best_new = new_brand_new_df.loc[best_new_idx].set_index('cluster_accession')

            pep_series = best_new['peptidoform'].astype(str)
            needs_charge = ~pep_series.str.contains('/', regex=False)
            pep_series = pep_series.where(
                ~needs_charge,
                pep_series + '/' + best_new['charge'].astype(int).astype(str))

            new_meta_df = pd.DataFrame({
                'cluster_id': best_new.index.astype(str),
                'species': species,
                'instrument': instrument,
                'charge': charge_int,
                'peptidoform': pep_series.values,
                'peptide_sequence': pep_series.str.split('/').str[0].str.replace(
                    r'\[.*?\]', '', regex=True).values,
                'consensus_mz_array': best_new['mz_array'].values,
                'consensus_intensity_array': best_new['intensity_array'].values,
                'consensus_method': method_type,
                'precursor_mz': best_new['pepmass'].astype(float).values,
                'member_count': 1,
                'project_count': 1,
                'best_pep': best_new['posterior_error_probability'].astype(float).values,
                'best_qvalue': best_new['global_qvalue'].astype(float).values,
                'purity': 1.0,
            })
            meta = pd.concat([existing_meta, new_meta_df], ignore_index=True)
    else:
        meta = existing_meta

    meta['cluster_id'] = meta['cluster_id'].astype(str)

    # Update stats from merged membership
    stats_map = cluster_stats.set_index('cluster_id')
    for col in ['member_count', 'project_count', 'best_pep', 'best_qvalue']:
        meta[col] = meta['cluster_id'].map(stats_map[col]).fillna(meta[col])
    meta['purity'] = meta['cluster_id'].map(purity_series).fillna(meta['purity'])

    schema = pa.schema([
        ('cluster_id', pa.string()),
        ('species', pa.string()),
        ('instrument', pa.string()),
        ('charge', pa.int8()),
        ('peptidoform', pa.string()),
        ('peptide_sequence', pa.string()),
        ('consensus_mz_array', pa.list_(pa.float32())),
        ('consensus_intensity_array', pa.list_(pa.float32())),
        ('consensus_method', pa.string()),
        ('precursor_mz', pa.float64()),
        ('member_count', pa.int32()),
        ('project_count', pa.int16()),
        ('best_pep', pa.float64()),
        ('best_qvalue', pa.float64()),
        ('purity', pa.float32()),
    ])

    # Ensure correct types
    meta['cluster_id'] = meta['cluster_id'].astype(str)
    meta['charge'] = meta['charge'].astype(int)
    meta['member_count'] = meta['member_count'].astype(int)
    meta['project_count'] = meta['project_count'].astype(int)
    meta['purity'] = meta['purity'].astype(float)

    table = pa.Table.from_pandas(meta[[f.name for f in schema]], schema=schema, preserve_index=False)
    pq_write.write_table(table, out_metadata, compression='zstd')

    click.echo(f"Incremental merge complete:")
    click.echo(f"  New clusters: {new_clusters_added}")
    click.echo(f"  Updated clusters: {existing_clusters_updated}")
    click.echo(f"  Total clusters: {len(meta)}")
    click.echo(f"  Metadata: {out_metadata}")
    click.echo(f"  Membership: {out_membership}")
