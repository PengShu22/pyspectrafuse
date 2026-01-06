from typing import Tuple
from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
import math
import copy
import logging

logger = logging.getLogger(__name__)


class BinningStrategy(ConsensusStrategy):
    """Binning strategy generates consensus spectrum by binning m/z values.
    
    This strategy bins m/z values across all spectra in a cluster and generates
    a new consensus spectrum by averaging intensities within each bin.
    """

    def __init__(self, min_mz: float = 100., max_mz: float = 2000., 
                 bin_size: float = 0.02, peak_quorum: float = 0.25, 
                 edge_case_threshold: float = 0.5):
        """Initialize BinningStrategy.
        
        Args:
            min_mz: Minimum m/z value to consider
            max_mz: Maximum m/z value to consider
            bin_size: Bin size in m/z units
            peak_quorum: Minimum fraction of spectra that must contain a peak
            edge_case_threshold: Threshold for handling edge cases at bin boundaries
        """
        self.min_mz: float = min_mz
        self.max_mz: float = max_mz
        self.bin_size: float = bin_size
        self.peak_quorum: float = peak_quorum
        self.edge_case_threshold: float = edge_case_threshold

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, 
                                       filter_metrics: str = 'global_qvalue') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Generate consensus spectrum using binning approach.
        
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

        res = merge_median_and_top.groupby('cluster_accession').apply(
            self.get_consensus_spectrum, include_groups=False)
        return res, single_spectrum_df

    def get_consensus_spectrum(self, single_group: pd.DataFrame) -> pd.Series:
        """Generate consensus spectrum for a single cluster group using binning.
        
        Args:
            single_group: DataFrame containing spectra from a single cluster
            
        Returns:
            Series representing the consensus spectrum
            
        Raises:
            ValueError: If single_group is empty or has no spectra
        """
        if single_group.empty:
            raise ValueError("single_group is empty")
            
        Nreps = single_group['Nreps'].iloc[0]
        posterior_error_probability = single_group['posterior_error_probability'].min()
        usi = ';'.join(single_group['usi'])

        # Create a dictionary for a cluster, key is usi, value is sus.MsmsSpectrum object
        cluster_spectra = single_group.set_index('usi')['ms2spectrumObj'].to_dict()
        
        if not cluster_spectra:
            raise ValueError("No spectra found in cluster")

        num_bins = math.ceil((self.max_mz - self.min_mz) / self.bin_size)
        mzs = np.zeros(num_bins, dtype=np.float32)
        intensities = np.zeros(num_bins, dtype=np.float32)
        n_peaks = np.zeros(num_bins, dtype=np.uint32)
        precursor_mzs, precursor_charges = [], []

        for spectrum in cluster_spectra.values():
            spectrum = copy.copy(spectrum).set_mz_range(self.min_mz, self.max_mz)
            bin_array = np.floor((spectrum.mz - self.min_mz) / self.bin_size).astype(np.uint32)
            # Ensure bin indices are within bounds
            valid_bins = bin_array < num_bins
            n_peaks[bin_array[valid_bins]] += 1
            intensities[bin_array[valid_bins]] += spectrum.intensity[valid_bins]
            mzs[bin_array[valid_bins]] += spectrum.mz[valid_bins]
            precursor_mzs.append(spectrum.precursor_mz)
            precursor_charges.append(spectrum.precursor_charge)

        # Handle edge cases where peaks are split on bin boundaries
        mz_temp = np.copy(mzs)
        mz_temp_mask = n_peaks > 0
        mz_temp[mz_temp_mask] /= n_peaks[mz_temp_mask]
        mz_delta = np.diff(mz_temp)
        mz_delta = np.append(mz_delta, 0)  # Add zero at end
        
        # Handle cases where deltas are smaller than the thresholded bin size
        mz_delta_small_index = np.nonzero(
            (mz_delta > 0) & (mz_delta < self.bin_size * self.edge_case_threshold))[0]
        if len(mz_delta_small_index) > 0:
            # Consolidate split peaks into one bin
            valid_indices = mz_delta_small_index[mz_delta_small_index + 1 < num_bins]
            mzs[valid_indices] += mzs[valid_indices + 1]
            mzs[valid_indices + 1] = 0
            intensities[valid_indices] += intensities[valid_indices + 1]
            intensities[valid_indices + 1] = 0
            n_peaks[valid_indices] += n_peaks[valid_indices + 1]
            n_peaks[valid_indices + 1] = 0

        # Determine how many peaks need to be present to keep a final peak
        peak_quorum_int = math.ceil(len(cluster_spectra) * self.peak_quorum)
        mask = n_peaks >= peak_quorum_int
        
        if not np.any(mask):
            raise ValueError("No peaks meet the quorum requirement")
        
        # Take the mean of all peaks per bin
        mzs = mzs[mask] / n_peaks[mask]
        intensities = intensities[mask] / n_peaks[mask]
        precursor_mz = np.mean(precursor_mzs)
        precursor_charge = precursor_charges[0]

        peptidoform = '; '.join(np.unique(single_group['peptidoform']))

        new_spectrum_index = ['pepmass', 'Nreps', 'posterior_error_probability', 'peptidoform',
                              'usi', 'charge', 'mz_array', 'intensity_array']
        return pd.Series([precursor_mz, Nreps, posterior_error_probability, peptidoform, 
                         usi, precursor_charge, mzs, intensities],
                         index=new_spectrum_index)
        Nreps = single_group['Nreps'].to_list()[0]
        posterior_error_probability = np.min(single_group['posterior_error_probability'])
        usi = ';'.join(single_group['usi'])

        # Create a dictionary for a cluster, key is usi, value is sus.MsmsSpectrum object
        cluster_spectra = single_group.set_index('usi')['ms2spectrumObj'].to_dict()
        spectra_keys = list(cluster_spectra.keys())

        num_bins = math.ceil((self.max_mz - self.min_mz) / self.bin_size)
        mzs = np.zeros(num_bins, dtype=np.float32)
        intensities = np.zeros(num_bins, dtype=np.float32)
        n_peaks = np.zeros(num_bins, dtype=np.uint32)
        precursor_mzs, precursor_charges = [], []

        for spectrum in cluster_spectra.values():
            spectrum = copy.copy(spectrum).set_mz_range(
                self.min_mz, self.max_mz)
            bin_array = np.floor((spectrum.mz - self.min_mz)
                                 / self.bin_size).astype(np.uint32)
            n_peaks[bin_array] += 1
            intensities[bin_array] += spectrum.intensity
            mzs[bin_array] += spectrum.mz
            precursor_mzs.append(spectrum.precursor_mz)
            precursor_charges.append(spectrum.precursor_charge)

        # # Verify that all precursor charges are the same.
        # if not all(charge == precursor_charges[0]
        #            for charge in precursor_charges):
        #     raise ValueError('Spectra in a cluster have different precursor '
        #                      'charges')
        # Try to handle the case where a single peak is split on a bin
        # boundary.
        mz_temp = np.copy(mzs)
        mz_temp_mask = n_peaks > 0
        mz_temp[mz_temp_mask] /= n_peaks[mz_temp_mask]
        # Subtract the mzs from their previous mz.
        mz_delta = np.diff(mz_temp)
        mz_delta[-1] = 0
        # Handle cases where the deltas are smaller than the thresholded bin
        # size.
        mz_delta_small_index = np.nonzero(
            (mz_delta > 0) &
            (mz_delta < self.bin_size * self.edge_case_threshold))[0]
        if len(mz_delta_small_index) > 0:
            # Consolidate all the split mzs, intensities, n_peaks into one bin.
            mzs[mz_delta_small_index] += mzs[mz_delta_small_index + 1]
            mzs[mz_delta_small_index + 1] = 0
            intensities[mz_delta_small_index] += \
                intensities[mz_delta_small_index + 1]
            intensities[mz_delta_small_index + 1] = 0
            n_peaks[mz_delta_small_index] += n_peaks[mz_delta_small_index + 1]
            n_peaks[mz_delta_small_index + 1] = 0

        # Determine how many peaks need to be present to keep a final peak.
        peak_quorum_int = math.ceil(len(cluster_spectra) * self.peak_quorum)
        mask = n_peaks >= peak_quorum_int
        # Take the mean of all peaks per bin.
        mzs = mzs[mask] / n_peaks[mask]
        intensities = intensities[mask] / n_peaks[mask]
        precursor_mz = np.mean(precursor_mzs)
        precursor_charge = precursor_charges[0]

        peptidoform = '; '.join(np.unique(single_group['peptidoform']))

        new_spectrum_index = ['pepmass', 'Nreps', 'posterior_error_probability', 'peptidoform',
                              'usi', 'charge', 'mz_array', 'intensity_array']
        # Return information for the newly generated consensus spectrum
        return pd.Series([precursor_mz, Nreps, posterior_error_probability, peptidoform, usi, precursor_charge, mzs, intensities],
                         index=new_spectrum_index)
