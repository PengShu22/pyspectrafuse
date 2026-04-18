from typing import Tuple, Callable
from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
import logging
import functools
import spectrum_utils.spectrum as sus
from pyspectrafuse.consensus_strategy.metrics import dot

logger = logging.getLogger(__name__)

# Clusters above this size use the faster centroid approximation instead
# of exact pairwise similarity. 50 covers >99% of real-world clusters.
_PAIRWISE_THRESHOLD = 50


class MostSimilarStrategy(ConsensusStrategy):
    """Most Similar strategy selects the spectrum most similar to all others in cluster.

    For small clusters (<=50 spectra): exact pairwise similarity via vectorized
    upper-triangle computation.
    For large clusters (>50 spectra): centroid approximation — compute mean
    spectrum, then find the real spectrum closest to it.
    """

    def __init__(self, sim: str = 'dot', fragment_mz_tolerance: float = 0.02):
        self.fragment_mz_tolerance: float = fragment_mz_tolerance
        self.sim: str = sim
        if self.sim == 'dot':
            self.compare_spectra: Callable = functools.partial(
                dot, fragment_mz_tolerance=fragment_mz_tolerance)
        else:
            raise ValueError(f"Unknown spectrum similarity method: {sim}. Only 'dot' is supported.")

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame,
                                       filter_metrics: str = 'global_qvalue') -> Tuple[pd.DataFrame, pd.DataFrame]:
        if cluster_df.empty:
            raise ValueError("cluster_df is empty")

        merge_median_and_top, single_spectrum_df = self.classify_cluster_group(cluster_df, filter_metrics)

        if merge_median_and_top.empty:
            logger.warning("No clusters with multiple spectra found")
            return merge_median_and_top, single_spectrum_df

        merge_median_and_top['ms2spectrumObj'] = merge_median_and_top.apply(
            lambda row: self.get_Ms2SpectrumObj(row), axis=1)

        res = merge_median_and_top.groupby('cluster_accession', group_keys=False).apply(
            self._select_representative)
        logger.info(f"Consensus spectrum shape: {res.shape}; Single spectrum shape: {single_spectrum_df.shape}")

        return res, single_spectrum_df

    def _select_representative(self, single_group: pd.DataFrame) -> pd.Series:
        """Select the most representative spectrum from a cluster group.

        Uses exact pairwise similarity for small clusters and centroid
        approximation for large ones.
        """
        if single_group.empty:
            raise ValueError("single_group is empty")

        n = len(single_group)
        if n == 1:
            return single_group.iloc[0]

        if n <= _PAIRWISE_THRESHOLD:
            return self._select_pairwise(single_group)
        else:
            return self._select_centroid(single_group)

    def _select_pairwise(self, single_group: pd.DataFrame) -> pd.Series:
        """Exact pairwise: compute upper triangle of similarity matrix."""
        cluster_spectra = single_group.set_index('usi')['ms2spectrumObj'].to_dict()
        spectra_keys = list(cluster_spectra.keys())
        n = len(spectra_keys)
        spectra_list = [cluster_spectra[key] for key in spectra_keys]

        # Pre-normalize intensities to avoid redundant normalization in dot()
        sim_matrix = np.zeros((n, n))

        # Compute upper triangle only
        for i in range(n):
            sim_matrix[i, i] = 1.0  # self-similarity
            for j in range(i + 1, n):
                similarity = self.compare_spectra(spectra_list[i], spectra_list[j])
                sim_matrix[i, j] = similarity
                sim_matrix[j, i] = similarity

        # Spectrum with highest total similarity to all others
        max_sim_index = sim_matrix.sum(axis=0).argmax()
        max_sim_usi = spectra_keys[max_sim_index]

        max_sim_row = single_group[single_group['usi'] == max_sim_usi]
        return max_sim_row.iloc[0]

    def _select_centroid(self, single_group: pd.DataFrame) -> pd.Series:
        """O(n) centroid approximation for large clusters.

        Computes a mean spectrum (centroid) from all members, then finds
        the real spectrum with the highest dot-product similarity to it.
        This is O(n) instead of O(n^2).
        """
        spectra_list = single_group['ms2spectrumObj'].tolist()
        usis = single_group['usi'].tolist()

        # Build centroid: bin all peaks, average intensities
        # Use the same binning as BinningStrategy but simpler
        all_mz = np.concatenate([s.mz for s in spectra_list])
        all_int = np.concatenate([s.intensity for s in spectra_list])

        if len(all_mz) == 0:
            return single_group.iloc[0]

        # Sort by m/z and cluster nearby peaks
        order = np.argsort(all_mz)
        all_mz = all_mz[order]
        all_int = all_int[order]

        # Merge peaks within fragment_mz_tolerance
        centroid_mz = []
        centroid_int = []
        i = 0
        while i < len(all_mz):
            j = i + 1
            while j < len(all_mz) and all_mz[j] - all_mz[i] <= self.fragment_mz_tolerance:
                j += 1
            # Average the cluster of peaks
            centroid_mz.append(np.mean(all_mz[i:j]))
            centroid_int.append(np.mean(all_int[i:j]))
            i = j

        centroid_mz = np.array(centroid_mz, dtype=np.float32)
        centroid_int = np.array(centroid_int, dtype=np.float32)

        # Create centroid as MsmsSpectrum
        precursor_mz = np.mean([s.precursor_mz for s in spectra_list])
        precursor_charge = spectra_list[0].precursor_charge
        centroid_spectrum = sus.MsmsSpectrum(
            'centroid', precursor_mz, precursor_charge, centroid_mz, centroid_int)

        # Find spectrum with highest similarity to centroid — O(n)
        best_idx = 0
        best_sim = -1.0
        for idx, spectrum in enumerate(spectra_list):
            sim = self.compare_spectra(centroid_spectrum, spectrum)
            if sim > best_sim:
                best_sim = sim
                best_idx = idx

        logger.debug(f"Centroid selection for cluster of {len(spectra_list)} spectra, "
                     f"best similarity: {best_sim:.4f}")

        return single_group.iloc[best_idx]
