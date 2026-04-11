"""Extract representative spectra from cluster_metadata.parquet to MGF.

Used for incremental clustering: representative spectra are re-clustered
alongside new data so that new spectra can merge into existing clusters.
"""
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

# MGF record template — avoids per-spectrum string concatenation overhead
_MGF_TEMPLATE = (
    "BEGIN IONS\n"
    "TITLE={title}\n"
    "PEPMASS={pepmass}\n"
    "CHARGE={charge}+\n"
    "{peaks}\n"
    "END IONS\n"
)


def extract_representatives(
    cluster_metadata_path: str,
    output_mgf_path: str,
    batch_size: int = 50_000,
    max_clusters: Optional[int] = None,
) -> int:
    """Write one representative spectrum per cluster to an MGF file.

    Each spectrum's TITLE is set to ``rep:{cluster_id}`` so that the
    representative can be traced back after re-clustering.

    Args:
        cluster_metadata_path: Path to an existing cluster_metadata.parquet.
        output_mgf_path: Destination MGF file path.
        batch_size: Arrow batch size for streaming reads.
        max_clusters: If set, only emit the first N clusters (useful for testing).

    Returns:
        Number of representative spectra written.
    """
    pf = pq.ParquetFile(cluster_metadata_path)
    columns = ['cluster_id', 'precursor_mz', 'charge',
               'consensus_mz_array', 'consensus_intensity_array']

    Path(output_mgf_path).parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with open(output_mgf_path, 'w') as fh:
        for batch in pf.iter_batches(batch_size=batch_size, columns=columns):
            df = batch.to_pandas()
            for _, row in df.iterrows():
                mz_arr = row['consensus_mz_array']
                int_arr = row['consensus_intensity_array']
                if mz_arr is None or int_arr is None:
                    continue
                if isinstance(mz_arr, np.ndarray):
                    mz_arr = mz_arr.tolist()
                if isinstance(int_arr, np.ndarray):
                    int_arr = int_arr.tolist()
                if len(mz_arr) == 0:
                    continue

                peaks = '\n'.join(
                    f'{mz} {intensity}'
                    for mz, intensity in zip(mz_arr, int_arr)
                )
                record = _MGF_TEMPLATE.format(
                    title=f"rep:{row['cluster_id']}",
                    pepmass=row['precursor_mz'],
                    charge=int(row['charge']),
                    peaks=peaks,
                )
                fh.write(record)
                written += 1

                if max_clusters is not None and written >= max_clusters:
                    logger.info(f"Reached max_clusters={max_clusters}, stopping")
                    return written

    logger.info(f"Wrote {written} representative spectra to {output_mgf_path}")
    return written
