"""Contains classes for filling trolley and fixed probe runs."""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import seaborn as sns

import gm2

class TrolleyRun(object):

    """Class for building DataFrames of Trolley runs."""

    def __init__(self, run):
        self.fp_run = gm2.FixedProbe([run])
        self.tr_run = gm2.Trolley([run])

        times, tr_freq_interp, tr_phi_interp, fp_freq_interp =\
            self.__time_interpolation()

        cols = ["tr_phi"]\
            + ["tr" + str(i) for i in np.arange(17)]\
            + ["fp" + str(i) for i in np.arange(378)]
        data = np.append(np.append(tr_phi_interp, tr_freq_interp, axis=1),
                         fp_freq_interp,
                         axis=1)
        self.interp_data = pd.DataFrame(data, index=times, columns=cols)

    def __time_interpolation(self):
        tr_time, tr_phi, tr_freq = self.tr_run.getBasics(mode_phi=2)
        tr_phi = np.rad2deg(np.unwrap(np.deg2rad(tr_phi)))
        _, fp_time, fp_freq = self.fp_run.getBasics()
        tr_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds
        fp_time /= 1.0e9

        tr_indices = np.mean(tr_freq, 1) > 0
        fp_indices = np.mean(fp_freq, 1) > 0

        times = np.arange(np.ceil(np.max([tr_time[tr_indices][0, 16],
                                          fp_time[fp_indices][0, 377]]
                                        )),
                          np.floor(np.min([fp_time[fp_indices][-1, 0],
                                           fp_time[fp_indices][-1, 0]]
                                         )) + 1, 1)

        tr_freq_interp = np.zeros([times.size, 17])
        tr_phi_interp = np.zeros([times.size, 1])
        fp_freq_interp = np.zeros([times.size, 378])

        for i in np.arange(17):
            tr_freq_interp[:, i] = interp1d(tr_time[tr_indices][:, i],
                                            tr_freq[tr_indices][:, i],
                                            kind='slinear')(times)
            tr_phi_interp[:, 0] += interp1d(tr_time[tr_indices][:, i],
                                            tr_phi[tr_indices][:, i],
                                            kind='slinear')(times)/17.
        for i in np.arange(378):
            fp_freq_interp[:, i] = interp1d(fp_time[fp_indices][:, i],
                                            fp_freq[fp_indices][:, i],
                                            kind='slinear')(times)

        return times, tr_freq_interp, tr_phi_interp, fp_freq_interp

    def probe_time_plot(self, probe='tr0'):
        """DOC STRING."""
        sns.set(style="darkgrid")

        sns.relplot(data=self.interp_data[probe], kind='line')
        fig = plt.gcf()
        axes = plt.gca()
        fig.set_size_inches(14, 8.66)
        xticks = axes.get_xticks()
        axes.set_xticklabels([pd.to_datetime(tm, unit='s').strftime('%Y-%m-%d\n %H:%M:%S')
                              for tm in xticks], rotation=0)
        return fig

    def probe_phi_plot(self, probe='tr0', wrap=True):
        """DOC STRING."""
        sns.set(style="darkgrid")

        pos_data = self.interp_data.loc[:, ['tr_phi', probe]].copy()
        if not wrap:
            pos_data['tr_phi'] = np.rad2deg(np.unwrap(np.deg2rad(pos_data['tr_phi'])))

        sns.relplot(x='tr_phi', y=probe, data=pos_data, kind='line')
        fig = plt.gcf()
        fig.set_size_inches(14, 8.66)

        return fig
