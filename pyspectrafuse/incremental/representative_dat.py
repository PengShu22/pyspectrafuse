"""Extract representative spectra from cluster_metadata.parquet to .dat format.

Used for incremental clustering: representative spectra are re-clustered
alongside new data so that new spectra can merge into existing clusters.

Replaces representative_mgf.py — writes MaRaCluster's binary .dat format
instead of MGF text, matching the full-mode pipeline path.
"""
import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pyarrow.parquet as pq

from pyspectrafuse.maracluster_dat import (
    _bin_spectrum, _pack_spectrum, _pack_scaninfo,
    PROTON_MASS,
)

logger = logging.getLogger(__name__)


def extract_representatives_dat(
    cluster_metadata_path: str,
    output_dir: str,
    charge_filter: Optional[int] = None,
    file_idx: int = 0,
    batch_size: int = 50_000,
    max_clusters: Optional[int] = None,
) -> Tuple[int, int, str, str, str]:
    """Write one representative spectrum per cluster to .dat binary format.

    Each spectrum's scan_title is set to ``rep:{cluster_id}`` so that the
    representative can be traced back after re-clustering.

    Args:
        cluster_metadata_path: Path to an existing cluster_metadata.parquet.
        output_dir: Directory to write .dat files.
        charge_filter: If set, only include clusters with this charge.
        file_idx: File index for the .dat struct (default 0).
        batch_size: Arrow batch size for streaming reads.
        max_clusters: If set, only emit the first N clusters (useful for testing).

    Returns:
        Tuple of (written, skipped, dat_path, scaninfo_path, titles_path).
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    stem = Path(cluster_metadata_path).stem
    if charge_filter is not None:
        stem = f"reps_{stem}_charge{charge_filter}"
    else:
        stem = f"reps_{stem}"

    dat_path = str(Path(output_dir) / f"{stem}.dat")
    scaninfo_path = str(Path(output_dir) / f"{stem}.scan_info.dat")
    titles_path = str(Path(output_dir) / f"{stem}.scan_titles.txt")

    pf = pq.ParquetFile(cluster_metadata_path)
    columns = ['cluster_id', 'precursor_mz', 'charge',
               'consensus_mz_array', 'consensus_intensity_array']

    written = 0
    skipped = 0

    with open(dat_path, 'wb') as f_dat, \
         open(scaninfo_path, 'wb') as f_scaninfo, \
         open(titles_path, 'w') as f_titles:

        for batch in pf.iter_batches(batch_size=batch_size, columns=columns):
            df = batch.to_pandas()

            if charge_filter is not None and 'charge' in df.columns:
                df = df[df['charge'] == charge_filter]
                if df.empty:
                    continue

            dat_buf = bytearray()
            scaninfo_buf = bytearray()
            title_lines = []

            for _, row in df.iterrows():
                mz_arr = row['consensus_mz_array']
                int_arr = row['consensus_intensity_array']
                if mz_arr is None or int_arr is None:
                    skipped += 1
                    continue

                mz_arr = np.asarray(mz_arr, dtype=np.float64)
                int_arr = np.asarray(int_arr, dtype=np.float64)
                if len(mz_arr) == 0:
                    skipped += 1
                    continue

                charge = int(row['charge'])
                if charge <= 0:
                    skipped += 1
                    continue

                prec_mz = float(row['precursor_mz'])
                prec_mass = prec_mz * charge - PROTON_MASS * (charge - 1)

                peak_bins = _bin_spectrum(mz_arr, int_arr, prec_mass)
                # Use min_peaks=1 for representatives — these are known-good
                # consensus spectra that should not be silently dropped
                if len(peak_bins) < 1:
                    skipped += 1
                    continue

                scannr = written
                cluster_id = str(row['cluster_id'])

                dat_buf.extend(_pack_spectrum(
                    file_idx, scannr, charge, prec_mz, 0.0, peak_bins))
                scaninfo_buf.extend(_pack_scaninfo(
                    file_idx, scannr, prec_mz, prec_mz))
                title_lines.append(
                    f"{file_idx}\t{scannr}\trep:{cluster_id}\n")

                written += 1
                if max_clusters is not None and written >= max_clusters:
                    break

            if dat_buf:
                f_dat.write(dat_buf)
                f_scaninfo.write(scaninfo_buf)
                f_titles.writelines(title_lines)

            if max_clusters is not None and written >= max_clusters:
                logger.info(f"Reached max_clusters={max_clusters}, stopping")
                break

    logger.info(f"Wrote {written} representative spectra to {dat_path} "
                f"({skipped} skipped)")
    return written, skipped, dat_path, scaninfo_path, titles_path
