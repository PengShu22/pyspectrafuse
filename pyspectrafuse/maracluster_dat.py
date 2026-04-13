"""Generate MaRaCluster-compatible .dat binary files directly from parquet.

Bypasses MGF entirely: reads QPX .psm.parquet → writes compact .dat binary
structs (100 bytes per spectrum). MaRaCluster can then skip its Step 1
(file conversion) and start directly from the pvalue computation step.

Size comparison for 1M PSMs:
  - MGF text:  ~3.4 GB
  - .dat binary: ~100 MB  (34x smaller)

The .dat format is a flat array of Spectrum structs:
  struct Spectrum {  // 100 bytes total
    uint32 fileIdx;         // 0-based file index
    uint32 scannr;          // scan number
    uint32 charge;          // precursor charge
    float  precMz;          // precursor m/z
    float  retentionTime;   // retention time
    int16  fragBins[40];    // top-40 binned peak indices, 0-padded
  };

Reference: github.com/statisticalbiotechnology/maracluster
  - src/Spectrum.h, src/BinSpectra.cpp, src/BinaryInterface.cpp
"""
import logging
import struct
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

# ── MaRaCluster constants ──
# Senko bin width (mass of one nucleon at ~1 Da resolution)
BIN_WIDTH = 1.000508
BIN_SHIFT = 0.32
MAX_SCORING_PEAKS = 40
MIN_SCORING_PEAKS = 15
PROTON_MASS = 1.00727646677

# Struct format: 2×uint32 (ScanId) + uint32 (charge) + 2×float32 + 40×int16
# Little-endian, total = 4+4+4+4+4+80 = 100 bytes
SPECTRUM_STRUCT = struct.Struct('<IIIff40h')
SCANINFO_STRUCT = struct.Struct('<IIff')

assert SPECTRUM_STRUCT.size == 100
assert SCANINFO_STRUCT.size == 16


def _bin_spectrum(mz_array: np.ndarray, intensity_array: np.ndarray,
                  prec_mass: float, n_peaks: int = MAX_SCORING_PEAKS) -> List[int]:
    """Replicate MaRaCluster's binBinaryTruncated algorithm.

    1. Sort peaks by intensity descending
    2. For each peak (highest intensity first):
       - Skip if mz >= precursor mass (neutral mass)
       - Compute bin index: floor(mz / 1.000508 + 0.32)
       - If bin not seen before, add to result
       - Stop at n_peaks unique bins
    3. Sort bins ascending

    Args:
        mz_array: Fragment m/z values
        intensity_array: Fragment intensities
        prec_mass: Neutral precursor mass (NOT m/z — charge-deconvolved)
        n_peaks: Maximum number of bins to keep (default 40)

    Returns:
        Sorted list of bin indices (may be shorter than n_peaks)
    """
    if len(mz_array) == 0:
        return []

    # Sort by intensity descending
    order = np.argsort(-intensity_array)

    seen = set()
    peak_bins = []
    for idx in order:
        mz = float(mz_array[idx])
        if mz >= prec_mass:
            continue
        b = int(mz / BIN_WIDTH + BIN_SHIFT)
        if b not in seen and b > 0:
            seen.add(b)
            peak_bins.append(b)
            if len(peak_bins) >= n_peaks:
                break

    peak_bins.sort()
    return peak_bins


def _pack_spectrum(file_idx: int, scannr: int, charge: int,
                   prec_mz: float, retention_time: float,
                   peak_bins: List[int]) -> bytes:
    """Pack one Spectrum struct into 100 bytes."""
    bins = list(peak_bins[:MAX_SCORING_PEAKS])
    bins.extend([0] * (MAX_SCORING_PEAKS - len(bins)))
    return SPECTRUM_STRUCT.pack(
        file_idx, scannr, charge,
        float(prec_mz), float(retention_time),
        *bins
    )


def _pack_scaninfo(file_idx: int, scannr: int,
                   min_prec_mz: float, max_prec_mz: float) -> bytes:
    """Pack one ScanInfo struct into 16 bytes."""
    return SCANINFO_STRUCT.pack(file_idx, scannr,
                                float(min_prec_mz), float(max_prec_mz))


def parquet_to_dat(
    parquet_path: str,
    output_dir: str,
    file_idx: int = 0,
    charge_filter: Optional[int] = None,
    batch_size: int = 50_000,
) -> Tuple[int, int, str]:
    """Convert a QPX .psm.parquet file to MaRaCluster .dat binary format.

    Writes three files:
      - {name}.dat          — array of 100-byte Spectrum structs
      - {name}.scan_info.dat — array of 16-byte ScanInfo structs
      - {name}.scan_titles.txt — TSV of (fileIdx, scannr, title)

    Args:
        parquet_path: Path to .psm.parquet file
        output_dir: Directory to write .dat files
        file_idx: 0-based file index for this file in the batch
        charge_filter: If set, only include spectra with this charge
        batch_size: Arrow batch size for streaming reads

    Returns:
        Tuple of (spectra_written, spectra_skipped, dat_path)
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    name = Path(parquet_path).stem
    if charge_filter is not None:
        name = f"{name}_charge{charge_filter}"
    dat_path = str(Path(output_dir) / f"{name}.dat")
    scaninfo_path = str(Path(output_dir) / f"{name}.scan_info.dat")
    titles_path = str(Path(output_dir) / f"{name}.scan_titles.txt")

    pf = pq.ParquetFile(parquet_path)
    columns = ['mz_array', 'intensity_array', 'charge',
               'observed_mz', 'calculated_mz', 'rt', 'scan',
               'run_file_name', 'peptidoform']

    # Read available columns (QPX format may use different names)
    schema_names = set(pf.schema_arrow.names)
    read_cols = [c for c in columns if c in schema_names]

    written = 0
    skipped = 0

    with open(dat_path, 'wb') as f_dat, \
         open(scaninfo_path, 'wb') as f_scaninfo, \
         open(titles_path, 'w') as f_titles:

        for batch in pf.iter_batches(batch_size=batch_size, columns=read_cols):
            df = batch.to_pandas()

            if charge_filter is not None and 'charge' in df.columns:
                df = df[df['charge'] == charge_filter]
                if df.empty:
                    continue

            for row_idx, (_, row) in enumerate(df.iterrows()):
                mz_arr = row.get('mz_array')
                int_arr = row.get('intensity_array')
                if mz_arr is None or int_arr is None:
                    skipped += 1
                    continue

                mz_arr = np.asarray(mz_arr, dtype=np.float64)
                int_arr = np.asarray(int_arr, dtype=np.float64)

                if len(mz_arr) == 0 or len(int_arr) == 0:
                    skipped += 1
                    continue

                charge = int(row.get('charge', 0))
                if charge <= 0:
                    skipped += 1
                    continue

                # Precursor m/z: prefer observed_mz, fall back to calculated_mz
                prec_mz = row.get('observed_mz')
                if prec_mz is None or np.isnan(prec_mz):
                    prec_mz = row.get('calculated_mz', 0.0)
                prec_mz = float(prec_mz)

                # Neutral mass for peak filtering
                prec_mass = prec_mz * charge - PROTON_MASS * (charge - 1)

                # Retention time — set to 0.0 to match MaRaCluster's MGF behavior
                # (MGF doesn't carry RT, so MaRaCluster always sees 0.0)
                rt = 0.0

                # Scan number — use sequential 0-based index, matching MaRaCluster's
                # MGF behavior (it assigns scan numbers 0, 1, 2, ... in file order)
                scannr = written

                # Bin the spectrum
                peak_bins = _bin_spectrum(mz_arr, int_arr, prec_mass)
                if len(peak_bins) < MIN_SCORING_PEAKS:
                    skipped += 1
                    continue

                # Write Spectrum struct
                f_dat.write(_pack_spectrum(
                    file_idx, scannr, charge, prec_mz, rt, peak_bins))

                # Write ScanInfo struct
                f_scaninfo.write(_pack_scaninfo(
                    file_idx, scannr, prec_mz, prec_mz))

                # Write title line — maps sequential index → original spectrum identity
                run = row.get('run_file_name', '')
                pep = row.get('peptidoform', '')
                orig_scan = row.get('scan')
                if orig_scan is not None:
                    if isinstance(orig_scan, (list, np.ndarray)):
                        orig_scan = int(orig_scan[0]) if len(orig_scan) > 0 else 0
                    else:
                        orig_scan = int(orig_scan)
                else:
                    orig_scan = 0
                title = f"id=mzspec::{run}:scan:{orig_scan}:{pep}/{charge}"
                f_titles.write(f"{file_idx}\t{scannr}\t{title}\n")

                written += 1

    dat_size = Path(dat_path).stat().st_size
    logger.info(f"Wrote {written} spectra to {dat_path} "
                f"({dat_size / 1e6:.1f} MB, {skipped} skipped)")

    return written, skipped, dat_path


def convert_dataset_to_dat(
    parquet_dir: str,
    output_dir: str,
    file_idx_start: int = 0,
    charge_filter: Optional[int] = None,
) -> Tuple[List[str], int]:
    """Convert all .psm.parquet files in a directory to .dat format.

    Args:
        parquet_dir: Directory containing .psm.parquet files
        output_dir: Directory to write .dat files
        file_idx_start: Starting file index
        charge_filter: If set, only include spectra with this charge

    Returns:
        Tuple of (list of .dat file paths, total spectra written)
    """
    from pyspectrafuse.commands.spectrum2msp import find_target_ext_files

    parquet_files = find_target_ext_files(parquet_dir, '.parquet')
    if not parquet_files:
        raise FileNotFoundError(f"No parquet files found in {parquet_dir}")

    dat_paths = []
    total_written = 0

    for i, pf in enumerate(parquet_files):
        written, skipped, dat_path = parquet_to_dat(
            pf, output_dir,
            file_idx=file_idx_start + i,
            charge_filter=charge_filter,
        )
        dat_paths.append(dat_path)
        total_written += written
        logger.info(f"  [{i+1}/{len(parquet_files)}] {Path(pf).name}: "
                     f"{written} written, {skipped} skipped")

    return dat_paths, total_written


def write_dat_file_list(dat_paths: List[str], output_path: str) -> str:
    """Write the file list that MaRaCluster's --datFNfile expects.

    Args:
        dat_paths: List of .dat file paths
        output_path: Where to write the file list

    Returns:
        output_path
    """
    with open(output_path, 'w') as f:
        for p in dat_paths:
            f.write(p + '\n')
    return output_path
