from typing import Tuple, Callable
from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
import logging
import functools
import spectrum_utils.spectrum as sus
from pyspectrafuse.consensus_strategy.metrics import dot

logger = logging.getLogger(__name__)


class MostSimilarStrategy(ConsensusStrategy):
    """Most Similar strategy selects the spectrum most similar to all others in cluster.
    
    This method calculates similarity distances between each spectrum in a cluster
    and all other spectra, sums the similarity distances, and selects the spectrum
    with the maximum sum as the representative spectrum.
    """

    def __init__(self, sim: str = 'dot', fragment_mz_tolerance: float = 0.02):
        """Initialize MostSimilarStrategy.
        
        Args:
            sim: Similarity measure method (currently only 'dot' supported)
            fragment_mz_tolerance: Fragment m/z tolerance for similarity calculation
            
        Raises:
            ValueError: If similarity method is not supported
        """
        self.fragment_mz_tolerance: float = fragment_mz_tolerance
        self.sim: str = sim
        if self.sim == 'dot':
            self.compare_spectra: Callable = functools.partial(
                dot, fragment_mz_tolerance=fragment_mz_tolerance)
        else:
            raise ValueError(f"Unknown spectrum similarity method: {sim}. Only 'dot' is supported.")

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame,
                                       filter_metrics: str = 'global_qvalue') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Generate consensus spectrum by selecting most similar spectrum per cluster.
        
        This method calculates similarity distances between each spectrum in a cluster
        and all other spectra, sums the similarity distances, and selects the spectrum
        with the maximum sum as the representative spectrum. The consensus spectrum
        comes from the spectrum itself (does not generate a new spectrum).
        
        Args:
            cluster_df: DataFrame containing cluster spectrum data
            filter_metrics: Metric column name for filtering (unused in this strategy)
            
        Returns:
            Tuple of (consensus_spectrum_df, single_spectrum_df)
            
        Raises:
            ValueError: If cluster_df is empty
        """
        if cluster_df.empty:
            raise ValueError("cluster_df is empty")
            
        merge_median_and_top, single_spectrum_df = self.classify_cluster_group(cluster_df, filter_metrics)
        
        if merge_median_and_top.empty:
            logger.warning("No clusters with multiple spectra found")
            return merge_median_and_top, single_spectrum_df

        merge_median_and_top['ms2spectrumObj'] = merge_median_and_top.apply(
            lambda row: self.get_Ms2SpectrumObj(row), axis=1)
        
        # Group by accession and apply selection function to each group
        res = merge_median_and_top.groupby('cluster_accession').apply(
            self._select_representative, include_groups=False)
        logger.info(f"Consensus spectrum shape: {res.shape}; Single spectrum shape: {single_spectrum_df.shape}")

        return res, single_spectrum_df

    def _select_representative(self, single_group: pd.DataFrame) -> pd.Series:
        """Select the most representative spectrum from a cluster group.
        
        Args:
            single_group: DataFrame containing spectra from a single cluster
            
        Returns:
            Series representing the most similar spectrum
            
        Raises:
            ValueError: If single_group is empty
        """
        if single_group.empty:
            raise ValueError("single_group is empty")
            
        # Create a dictionary for a cluster, key is usi, value is sus.MsmsSpectrum object
        cluster_spectra = single_group.set_index('usi')['ms2spectrumObj'].to_dict()
        spectra_keys = list(cluster_spectra.keys())
        
        if len(spectra_keys) == 1:
            return single_group.iloc[0]

        # Build similarity matrix - optimized for symmetric matrix
        # Only calculate upper triangle and mirror to lower triangle
        n = len(spectra_keys)
        sim_matrix = np.zeros((n, n))
        spectra_list = [cluster_spectra[key] for key in spectra_keys]
        
        # Calculate upper triangle including diagonal
        for i in range(n):
            for j in range(i, n):
                similarity = self.compare_spectra(spectra_list[i], spectra_list[j])
                sim_matrix[i, j] = similarity
                if i != j:  # Mirror to lower triangle (skip diagonal)
                    sim_matrix[j, i] = similarity
        
        # Find the spectrum with the maximum similarity to all other spectra
        # Use sum along axis=0 to get total similarity for each spectrum
        max_sim_index = sim_matrix.sum(axis=0).argmax()
        max_sim_usi = spectra_keys[max_sim_index]

        # Find and return the row corresponding to the most similar spectrum
        max_sim_row = single_group[single_group['usi'] == max_sim_usi]
        return max_sim_row.iloc[0]
