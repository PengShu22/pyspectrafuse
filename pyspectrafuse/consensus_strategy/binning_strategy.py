from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
import math
import copy


class BinningStrategy(ConsensusStrategy):
    def __init__(self, min_mz=100., max_mz=2000., bin_size=0.02, peak_quorum=0.25, edge_case_threshold=0.5):
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.bin_size = bin_size
        self.peak_quorum = peak_quorum
        self.edge_case_threshold = edge_case_threshold

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, filter_metrics='global_qvalue'):

        merge_median_and_top, single_spectrum_df = self.classify_cluster_group(cluster_df, filter_metrics)
        merge_median_and_top['ms2spectrumObj'] = merge_median_and_top.apply(lambda row:
                                                                            self.get_Ms2SpectrumObj(row), axis=1)

        res = merge_median_and_top.groupby('cluster_accession').apply(self.get_consensus_spectrum)
        return res, single_spectrum_df

    def get_consensus_spectrum(self, single_group):
        Nreps = single_group['Nreps'].to_list()[0]
        posterior_error_probability = np.min(single_group['posterior_error_probability'])
        usi = ';'.join(single_group['usi'])

        # 内部应该为一个簇生成一个字典，key为usi， value为sus.MsmsSpectrum对象
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

        new_spectrum_index = ['pepmass', 'Nreps', 'posterior_error_probability',
                              'usi', 'charge', 'mz_array', 'intensity_array']
        # 返回生成的新的共识谱的信息
        return pd.Series([precursor_mz, Nreps, posterior_error_probability, usi, precursor_charge, mzs, intensities],
                         index=new_spectrum_index)  # 需要加usi, 然后comment用
