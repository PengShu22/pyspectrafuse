"""QPX metadata reader — reads .run.parquet and .sample.parquet.

Replaces SdrfUtil for runtime metadata lookups. Returns the same
{run_file_name: [species, instrument]} dict interface so replacement
is seamless across the pipeline.
"""

import glob
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import pyarrow.parquet as pq

logger = logging.getLogger(__name__)


def find_qpx_metadata_files(folder: str):
    """Find .run.parquet and .sample.parquet files in a directory.

    Uses glob.glob with recursive=True which follows symlinks (unlike
    Path.rglob which does not). This is needed for Nextflow staging dirs.

    Returns:
        (run_paths, sample_paths): Lists of paths found.
    """
    run_paths = sorted(glob.glob(str(Path(folder) / "**/*.run.parquet"), recursive=True))
    sample_paths = sorted(glob.glob(str(Path(folder) / "**/*.sample.parquet"), recursive=True))
    return run_paths, sample_paths


def read_run_parquet(run_path: str) -> dict:
    """Read a .run.parquet file and return {run_file_name: instrument} dict."""
    table = pq.read_table(run_path, columns=["run_file_name", "instrument"])
    df = table.to_pandas()
    return dict(zip(df["run_file_name"], df["instrument"].fillna("Unknown")))


def read_sample_parquet(sample_path: str) -> dict:
    """Read a .sample.parquet file and return {sample_accession: organism} dict."""
    table = pq.read_table(sample_path, columns=["sample_accession", "organism"])
    df = table.to_pandas()
    return dict(zip(df["sample_accession"], df["organism"].fillna("Unknown")))


def read_run_sample_link(run_path: str) -> dict:
    """Read .run.parquet and return {run_file_name: sample_accession} mapping.

    For multiplexed runs with multiple samples, returns the first sample.
    """
    table = pq.read_table(run_path, columns=["run_file_name", "samples"])
    result = {}
    for i in range(table.num_rows):
        run_file = table.column("run_file_name")[i].as_py()
        samples = table.column("samples")[i].as_py()
        if samples:
            result[run_file] = samples[0]["sample_accession"]
        else:
            result[run_file] = None
    return result


def get_metadata_dict_from_qpx(run_path: str, sample_path: str,
                                skip_instrument: bool = False) -> Dict[str, List[str]]:
    """Build {run_file_name: [species, instrument]} dict from QPX parquets.

    Drop-in replacement for SdrfUtil.get_metadata_dict_from_sdrf().

    Args:
        run_path: Path to .run.parquet file.
        sample_path: Path to .sample.parquet file.
        skip_instrument: If True, use 'all_instruments' for all runs.

    Returns:
        Dict mapping run_file_name to [organism, instrument].
    """
    # Read run -> instrument mapping
    run_instrument = read_run_parquet(run_path)

    # Read sample -> organism mapping
    sample_organism = read_sample_parquet(sample_path)

    # Read run -> sample link
    run_sample = read_run_sample_link(run_path)

    # Build the combined dict
    result = {}
    for run_file, sample_acc in run_sample.items():
        organism = sample_organism.get(sample_acc, "Unknown") if sample_acc else "Unknown"
        instrument = "all_instruments" if skip_instrument else run_instrument.get(run_file, "Unknown")
        result[run_file] = [organism, instrument]

    return result


def get_metadata_dict(folder: str, skip_instrument: bool = False,
                      run_paths: Optional[List[str]] = None,
                      sample_paths: Optional[List[str]] = None) -> Dict[str, List[str]]:
    """Get {run_file_name: [species, instrument]} from QPX metadata in a folder.

    Merges metadata from all .run.parquet and .sample.parquet files found.

    Args:
        folder: Directory containing QPX parquet files.
        skip_instrument: If True, use 'all_instruments' for all runs.
        run_paths: Optional explicit list of run parquet paths (overrides folder scan).
        sample_paths: Optional explicit list of sample parquet paths (overrides folder scan).

    Returns:
        Dict mapping run_file_name to [organism, instrument].
    """
    if run_paths is None or sample_paths is None:
        found_runs, found_samples = find_qpx_metadata_files(folder)
        if run_paths is None:
            run_paths = found_runs
        if sample_paths is None:
            sample_paths = found_samples

    if not run_paths:
        raise FileNotFoundError(f"No .run.parquet files found in {folder}")
    if not sample_paths:
        raise FileNotFoundError(f"No .sample.parquet files found in {folder}")

    # Merge all run and sample files
    merged = {}
    for run_path, sample_path in zip(run_paths, sample_paths):
        d = get_metadata_dict_from_qpx(run_path, sample_path,
                                         skip_instrument=skip_instrument)
        merged.update(d)

    # Handle case where more run files than sample files (or vice versa)
    if len(run_paths) != len(sample_paths):
        # Build complete sample dict from all sample files
        all_samples = {}
        for sp in sample_paths:
            all_samples.update(read_sample_parquet(sp))

        # Rebuild with all run files against merged sample dict
        merged = {}
        for rp in run_paths:
            run_instrument = read_run_parquet(rp)
            run_sample = read_run_sample_link(rp)
            for run_file, sample_acc in run_sample.items():
                organism = all_samples.get(sample_acc, "Unknown") if sample_acc else "Unknown"
                instrument = "all_instruments" if skip_instrument else run_instrument.get(run_file, "Unknown")
                merged[run_file] = [organism, instrument]

    logger.info(f"Loaded QPX metadata for {len(merged)} runs")
    return merged
