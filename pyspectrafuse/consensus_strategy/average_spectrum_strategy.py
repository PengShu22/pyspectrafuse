from pyteomics import mass
from pyspectrafuse.consensus_strategy.consensus_strategy_base import ConsensusStrategy
import pandas as pd
import numpy as np
from pyspectrafuse.mgf_convert.parquet2mgf import Parquet2Mgf
from pyspectrafuse.common.msp_utils import MspUtil


class AverageSpectrumStrategy(ConsensusStrategy):
    H = mass.nist_mass['H+'][0][0]

    def __init__(self, DIFF_THRESH=0.01, DYN_RANGE=1000, MIN_FRACTION=0.5, pepmass='lower_median', msms_avg='weighted'):
        self.DIFF_THRESH = DIFF_THRESH
        self.DYN_RANGE = DYN_RANGE
        self.MIN_FRACTION = MIN_FRACTION
        self.pepmass = pepmass
        self.msms_avg = msms_avg

        msms_avg_list = ['naive', 'weighted']
        pepmass_lst = ['naive_average', 'neutral_average', 'lower_median']

        self.params_dict = {
            "mode": ['encoded_clusters'],
            'dyn_range': DYN_RANGE,
            'min_fraction': MIN_FRACTION,
            'mz_accuracy': DIFF_THRESH,
            'pepmass': self.pepmass,  # lower_median
            'msms_avg': self.msms_avg  # weighted
        }
        self.get_pepmass = {'naive_average': self.naive_average_mass_and_charge,
                            'neutral_average': self.neutral_average_mass_and_charge,
                            'lower_median': self.lower_median_mass}[pepmass]
        self.kwargs = {'mz_accuracy': self.params_dict['mz_accuracy'], 'dyn_range': self.params_dict['dyn_range'],
                       'min_fraction': self.params_dict['min_fraction'], 'msms_avg': self.params_dict['msms_avg']}

    def consensus_spectrum_aggregation(self, cluster_df: pd.DataFrame, filter_metrics='global_qvalue'):
        merge_median_and_top, single_spectrum_df = self.classify_cluster_group(cluster_df, filter_metrics)
        merge_median_and_top['ms2spectrumDict'] = merge_median_and_top.apply(lambda row: self.get_Ms2SpectrumDict(row),
                                                                             axis=1)
        res = merge_median_and_top.groupby('cluster_accession').apply(self.get_average_spectrum)

        return res, single_spectrum_df

    def average_spectrum(self, spectra, pepmass='', charge='', Nreps='', pep='', usi=''):
        '''
        Produces an average spectrum for a cluster.

        Parameters
        ----------
        spectra : Iterable
            An iterable of spectra, each spectrum is expected to be in pyteomics format.
        title : str, optional
            Title of the output spectrum
        mz_accuracy : float, keyword only, optional
            m/z accuracy used for MS/MS clustering. Default is `DIFF_THRESH`
        dyn_range : float, keyword only, optional
            Dynamic range to apply to output (peaks less than max_intensity / dyn_range
            are discarded)
        msms_avg : {'naive', 'weighted'}
            Method for calculation of MS/MS peak m/z values in representative spectrum.
            Naive: simple average of peak m/z within MS/MS-level cluster.
            Weighted: weighted average with MS/MS peak intensities as weights.
        min_fraction, float, keyword only, optional
            Minimum fraction of cluster spectra need to contain the peak.
        Nreps： cluster

        Returns
        -------
        out : average spectrum in pyteomics format


                生成一个聚类中的平均频谱。

            参数
            ----------
            spectra:可迭代的
            光谱的可迭代对象，每个光谱期望为pyteomics格式。
            Title: str，可选
            输出频谱的标题
            Mz_accuracy:浮点数，仅限关键字，可选
            m/z精度用于MS/MS聚类。默认值是`DIFF_THRESH`
            Dyn_range:浮点数，仅限关键字，可选
            应用于输出的动态范围(峰值小于max_intensity / dyn_range
            被丢弃)
            Msms_avg: {'naive'， 'weighted'}
            代表光谱中MS/MS峰值m/z值的计算方法。
            Naive: MS/MS级别集群内峰值m/z的简单平均值。
            加权:以MS/MS峰值强度为权重的加权平均。
            Min_fraction，浮点数，仅限关键字，可选
            簇谱的最小部分需要包含峰。

            返回
            -------
            输出:pyteomics格式的平均光谱
            :param usi:
            :param pep:
            :param Nreps:
        '''
        mz_accuracy = self.kwargs.get('mz_accuracy', self.DIFF_THRESH)
        dyn_range = self.kwargs.get('dyn_range', self.DYN_RANGE)
        min_fraction = self.kwargs.get('min_fraction', self.MIN_FRACTION)
        msms_avg = self.kwargs.get('msms_avg')

        mz_arrays, int_arrays = [], []
        n = 0  # number of spectra
        for s in spectra:
            mz_arrays.append(s['m/z array'])
            int_arrays.append(s['intensity array'])
            n += 1
        if n > 1:
            mz_array_all = np.concatenate(mz_arrays)
            intensity_array_all = np.concatenate(int_arrays)

            idx = np.argsort(mz_array_all)
            mz_array_all = mz_array_all[idx]
            intensity_array_all = intensity_array_all[idx]
            diff_array = np.diff(mz_array_all)

            new_mz_array = []
            new_intensity_array = []

            ind_list = list(np.where(diff_array >= mz_accuracy)[0] + 1)

            i_prev = ind_list[0]

            mz_array_sum = np.cumsum(mz_array_all)
            intensity_array_sum = np.cumsum(intensity_array_all)
            if msms_avg == 'weighted':
                mult_mz_intensity_array_sum = np.cumsum(mz_array_all * intensity_array_all)

            min_l = min_fraction * n
            if i_prev >= min_l:
                I_sum = intensity_array_sum[i_prev - 1]
                if msms_avg == 'naive':
                    new_mz_array.append(mz_array_sum[i_prev - 1] / i_prev)
                elif msms_avg == 'weighted':
                    mzI_sum = mult_mz_intensity_array_sum[i_prev - 1]
                    new_mz_array.append(mzI_sum / I_sum)
                new_intensity_array.append(I_sum / n)

            for i in ind_list[1:-1]:
                if i - i_prev >= min_l:
                    I_sum = intensity_array_sum[i - 1] - intensity_array_sum[i_prev - 1]
                    if msms_avg == 'naive':
                        mz_sum = mz_array_sum[i - 1] - mz_array_sum[i_prev - 1]
                        new_mz_array.append(mz_sum / (i - i_prev))
                    elif msms_avg == 'weighted':
                        mzI_sum = mult_mz_intensity_array_sum[i - 1] - mult_mz_intensity_array_sum[i_prev - 1]
                        new_mz_array.append(mzI_sum / I_sum)
                    new_intensity_array.append(I_sum / n)
                i_prev = i

            if (len(mz_array_sum) - i_prev) >= min_l:
                I_sum = intensity_array_sum[-1] - intensity_array_sum[i_prev - 1]
                if msms_avg == 'naive':
                    mz_sum = mz_array_sum[-1] - mz_array_sum[i_prev - 1]
                    new_mz_array.append(mz_sum / (len(mz_array_sum) - i_prev))
                elif msms_avg == 'weighted':
                    mzI_sum = mult_mz_intensity_array_sum[-1] - mult_mz_intensity_array_sum[i_prev - 1]
                    new_mz_array.append(mzI_sum / I_sum)
                new_intensity_array.append(I_sum / n)
        else:
            new_mz_array = mz_arrays[0]
            new_intensity_array = int_arrays[0]

        new_mz_array = np.array(new_mz_array)
        new_intensity_array = np.array(new_intensity_array)

        min_i = new_intensity_array.max() / dyn_range
        idx = new_intensity_array >= min_i
        new_intensity_array = new_intensity_array[idx]
        new_mz_array = new_mz_array[idx]

        new_spectrum_index = ['pepmass', 'Nreps', 'posterior_error_probability',
                              'usi', 'charge', 'mz_array', 'intensity_array']
        # 返回生成的新的共识谱的信息
        return pd.Series([pepmass, Nreps, pep, usi, charge, new_mz_array, new_intensity_array],
                         index=new_spectrum_index)

    @staticmethod
    def _lower_median_mass_index(masses):
        i = np.argsort(masses)
        k = (len(masses) - 1) // 2
        idx = i[k]
        return idx, masses[idx]

    def lower_median_mass(self, spectra):
        masses, charges = self._neutral_masses(spectra)
        i, m = self._lower_median_mass_index(masses)
        z = charges[i]
        return (m + z * self.H) / z, z

    def lower_median_mass_rt(self, spectra):
        masses, charges = self._neutral_masses(spectra)
        rts = [s['params']['rtinseconds'] for s in spectra]
        i, m = self._lower_median_mass_index(masses)
        return rts[i]

    # def get_cluster_id(title):
    #     return title.split(';', 1)[0]

    @staticmethod
    def naive_average_mass_and_charge(spectra):
        mzs = [s['params']['pepmass'][0] for s in spectra]
        charges = {tuple(s['params']['charge']) for s in spectra}
        if len(charges) > 1:
            raise ValueError('There are different charge states in the cluster. Cannot average precursor m/z.')
        return sum(mzs) / len(mzs), charges.pop()[0]

    def _neutral_masses(self, spectra):
        mzs = [s['params']['pepmass'][0] for s in spectra]
        charges = [s['params']['charge'][0] for s in spectra if len(s['params']['charge']) == 1]
        masses = [(m * c - c * self.H) for m, c in zip(mzs, charges)]
        return masses, charges

    def neutral_average_mass_and_charge(self, spectra):
        masses, charges = self._neutral_masses(spectra)
        z = int(round(sum(charges) / len(charges)))
        avg_mass = sum(masses) / len(masses)
        return (avg_mass + z * self.H) / z, z

    @staticmethod
    def median_rt(spectra):
        rts = [s['params']['rtinseconds'] for s in spectra]
        return np.median(rts)

    def get_average_spectrum(self, single_group):
        Nreps = single_group['Nreps'].to_list()[0]
        posterior_error_probability = np.min(single_group['posterior_error_probability'])
        usi = ';'.join(single_group['usi'])

        spectra = single_group['ms2spectrumDict'].to_list()
        single_spectra_params = spectra[0]['params']
        pepmass_spectrum, charge = self.get_pepmass(spectra)
        return self.average_spectrum(spectra=spectra,
                                     pepmass=pepmass_spectrum,
                                     charge=charge,
                                     Nreps=Nreps,
                                     pep=posterior_error_probability,
                                     usi=usi)

