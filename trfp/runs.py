"""Contains classes for filling trolley and fixed probe runs."""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import seaborn as sns

import gm2
import trfp

from datetime import datetime
import pytz


class TrolleyRun(object):

    """Class for building DataFrames of Trolley runs."""

    def __init__(self, run):
        self.run = run
        self.fp_run = gm2.FixedProbe([run])
        self.tr_run = gm2.Trolley([run])

        times, tr_freq_interp, tr_phi_interp, fp_freq_interp =\
            self.__time_interpolation()
        
        # apply plunging probe calibrations to tr_freq_interp
        for ii in np.arange(17):
            tr_freq_interp[:,ii] = tr_freq_interp[:,ii] + trfp.PLUNGING_PROBE_CALIBRATIONS[ii]

        cols = ["tr_phi"]\
            + ["tr" + str(i) for i in np.arange(17)]\
            + ["fp" + str(i) for i in np.arange(378)]
        data = np.append(np.append(tr_phi_interp, tr_freq_interp, axis=1),
                         fp_freq_interp,
                         axis=1)
        self.interp_df = pd.DataFrame(data, index=times, columns=cols)
        
        self.moment_df = self.__moment_dataframe()

    def __time_interpolation(self):
        tr_time, tr_phi, tr_freq = self.tr_run.getBasics(mode_phi=2)
        tr_phi = np.rad2deg(np.unwrap(np.deg2rad(tr_phi)))
        _, fp_time, fp_freq = self.fp_run.getBasics()
        tr_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds
        fp_time /= 1.0e9

        tr_indices = np.mean(tr_freq, 1) > 0
        fp_indices = np.mean(fp_freq, 1) > 0

        times = np.arange(np.ceil(np.max([tr_time[tr_indices][0, 16],
                                          fp_time[fp_indices][0, 377]])),
                          np.floor(np.min([fp_time[fp_indices][-1, 0],
                                           fp_time[fp_indices][-1, 0]])) + 1,
                          1)

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

    def __moment_dataframe(self):
        """Builds dataframe of field moments."""

        # copy the interpolation data for the phi
        moment_df = self.interp_df['tr_phi'].copy()

        # create the 17 trolley moments at each point in time
        print 'Calculating trolley moments.',
        temp_df = self.interp_df[['tr' + str(tr)
                                    for tr in np.arange(17)]].copy()
        tr_st_m = np.zeros([temp_df.shape[0], 17])

        for index in np.arange(temp_df.shape[0]):
            tr_st_m[index, :] = np.dot(trfp.THETA_TR, temp_df.iloc[index])

        moment_df = pd.concat([moment_df,
                               pd.DataFrame(tr_st_m, index=temp_df.index,
                                            columns=['tr,m'+str(fp_m)
                                                     for fp_m in np.arange(1, 18)])],
                              axis=1, join='inner')

        # create the 72*6 fixed probe moments
        for station in np.arange(72):
            print '\rCalculating station ' + str(station) + ' moments.',
            temp_df = self.interp_df[['fp'+str(fp)
                                        for fp
                                        in trfp.STATION_PROBE_ID[station]]].copy()
            fp_st_m = np.zeros([temp_df.shape[0], 6])
            for index in np.arange(temp_df.shape[0]):
                if temp_df.shape[1] == 6:
                    fp_st_m[index, :] = np.dot(trfp.THETA_FP_6, temp_df.iloc[index])
                elif temp_df.shape[1] == 4:
                    fp_st_m[index, 0:4] = np.dot(trfp.THETA_FP_4, temp_df.iloc[index])
                    fp_st_m[index, 4:6] = [np.nan, np.nan]
            moment_df = pd.concat([moment_df,
                                   pd.DataFrame(fp_st_m, index=temp_df.index,
                                                columns=['st' + str(station) + ',m' + str(fp_m)
                                                         for fp_m in np.arange(1, 7)])],
                                  axis=1, join='inner')
        print '\rFinished calculating all moments for ' + str(temp_df.shape[0]) + ' events.'

        return moment_df
    
    def moment_time_plot(self, moment='tr,m1', save=False, ppm=True):
        """DOC STRING."""
        sns.set(style="darkgrid")
        
        if ppm == True:
            fg = sns.relplot(data=self.moment_df[moment]/61.79, kind='line')
        else:
            fg = sns.relplot(data=self.moment_df[moment], kind='line')
        fig = plt.gcf()
        for ax in fg.axes.flat:
            xticks = ax.get_xticks()
            ax.set_xticklabels([pd.to_datetime(tm, unit='s').tz_localize('UTC').tz_convert('US/Central').strftime('%Y-%m-%d\n %H:%M:%S %Z')
                                  for tm in xticks], rotation=0)
        fig.set_size_inches(6.5, 4)

        return fig
    
    def moment_correlation_plot(self, x_moment='tr,m1', y_moment='tr,m2', save=False, ppm=True):
        """DOC STRING."""
        sns.set(style="darkgrid")
        if ppm ==True:
            fg = sns.regplot(x_moment, y_moment, data=self.moment_df/61.79)
        else:
            fg = sns.regplot(x_moment, y_moment, data=self.moment_df)
        fig = plt.gcf()
        fig.set_size_inches(6.5, 4)

        return fig    

    def probe_time_plot(self, probe='tr0'):
        """DOC STRING."""
        sns.set(style="darkgrid")

        ax = sns.relplot(data=self.interp_df[probe], kind='line')
        fig = plt.gcf()
        fig.set_size_inches(14, 8.66)
        xticks = ax.get_xticks()
        ax.set_xticklabels([pd.to_datetime(tm, unit='s').tz_localize('UTC').tz_convert('US/Central').strftime('%Y-%m-%d\n %H:%M:%S %Z')
                              for tm in xticks], rotation=0)
        return fig

    def probe_phi_plot(self, probe='tr0', wrap=True):
        """DOC STRING."""
        sns.set(style="darkgrid")

        pos_data = self.interp_df.loc[:, ['tr_phi', probe]].copy()
        if not wrap:
            pos_data['tr_phi'] = np.rad2deg(np.unwrap(np.deg2rad(pos_data['tr_phi'])))

        sns.relplot(x='tr_phi', y=probe, data=pos_data, kind='line')
        fig = plt.gcf()
        fig.set_size_inches(14, 8.66)

        return fig

    
 
class FixedProbeRun(object):

    """Class for building DataFrames of fixed frobe runs."""

    def __init__(self, run):
        self.run = run
        self.fp_run = gm2.FixedProbe([run])
        # self.tr_run = gm2.Trolley([run])

        times, fp_freq_interp =\
            self.__time_interpolation()

        cols = ["fp" + str(i) for i in np.arange(378)]
        data = fp_freq_interp
        self.interp_df = pd.DataFrame(data, index=times, columns=cols)
        
        self.moment_df = self.__moment_dataframe()

    def __time_interpolation(self):
        _, fp_time, fp_freq = self.fp_run.getBasics()
        fp_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds

        fp_indices = np.mean(fp_freq, 1) > 0

        times = np.arange(np.ceil(fp_time[fp_indices][0, 377]),
                          np.floor(fp_time[fp_indices][-1, 0]) + 1,
                          2)

        fp_freq_interp = np.zeros([times.size, 378])

        for i in np.arange(378):
            fp_freq_interp[:, i] = interp1d(fp_time[fp_indices][:, i],
                                            fp_freq[fp_indices][:, i],
                                            kind='slinear')(times)

        return times, fp_freq_interp

    def __moment_dataframe(self):
        """Builds dataframe of field moments."""

        # copy the interpolation data for the index
        station = 0
        print '\rCalculating station ' + str(station) + ' moments.',
        temp_df = self.interp_df[['fp'+str(fp)
                                    for fp
                                    in trfp.STATION_PROBE_ID[station]]].copy()
        fp_st_m = np.zeros([temp_df.shape[0], 6])
        for index in np.arange(temp_df.shape[0]):
            if temp_df.shape[1] == 6:
                fp_st_m[index, :] = np.dot(trfp.THETA_FP_6, temp_df.iloc[index])
            elif temp_df.shape[1] == 4:
                fp_st_m[index, 0:4] = np.dot(trfp.THETA_FP_4, temp_df.iloc[index])
                fp_st_m[index, 4:6] = [np.nan, np.nan]
        moment_df = pd.DataFrame(fp_st_m, index=temp_df.index,
                                 columns=['st' + str(station) + ',m' + str(fp_m)
                                          for fp_m in np.arange(1, 7)])

        # create the 72*6 fixed probe moments
        for station in np.arange(1,72):
            print '\rCalculating station ' + str(station) + ' moments.',
            temp_df = self.interp_df[['fp'+str(fp)
                                        for fp
                                        in trfp.STATION_PROBE_ID[station]]].copy()
            fp_st_m = np.zeros([temp_df.shape[0], 6])
            for index in np.arange(temp_df.shape[0]):
                if temp_df.shape[1] == 6:
                    fp_st_m[index, :] = np.dot(trfp.THETA_FP_6, temp_df.iloc[index])
                elif temp_df.shape[1] == 4:
                    fp_st_m[index, 0:4] = np.dot(trfp.THETA_FP_4, temp_df.iloc[index])
                    fp_st_m[index, 4:6] = [np.nan, np.nan]
            moment_df = pd.concat([moment_df,
                                   pd.DataFrame(fp_st_m, index=temp_df.index,
                                                columns=['st' + str(station) + ',m' + str(fp_m)
                                                         for fp_m in np.arange(1, 7)])],
                                  axis=1, join='inner')
        print '\rFinished calculating all moments for ' + str(temp_df.shape[0]) + ' events.'

        return moment_df
    
    def moment_time_plot(self, moment='st30,m1', save=False, ppm=True):
        """DOC STRING."""
        sns.set(style="darkgrid")
        
        if ppm == True:
            fg = sns.relplot(data=self.moment_df[moment]/61.79, kind='line')
        else:
            fg = sns.relplot(data=self.moment_df[moment], kind='line')
        fig = plt.gcf()
        for ax in fg.axes.flat:
            xticks = ax.get_xticks()
            ax.set_xticklabels([pd.to_datetime(tm, unit='s').tz_localize('UTC').tz_convert('US/Central').strftime('%Y-%m-%d\n %H:%M:%S %Z')
                                  for tm in xticks], rotation=0)
        fig.set_size_inches(6.5, 4)

        return fig
    
    def moment_correlation_plot(self, x_moment='st30,m1', y_moment='st66,m1', save=False, ppm=True):
        """DOC STRING."""
        sns.set(style="darkgrid")
        if ppm ==True:
            fg = sns.regplot(x_moment, y_moment, data=self.moment_df/61.79)
        else:
            fg = sns.regplot(x_moment, y_moment, data=self.moment_df)
        fig = plt.gcf()
        fig.set_size_inches(6.5, 4)

        return fig    

    def probe_time_plot(self, probe='tr0'):
        """DOC STRING."""
        sns.set(style="darkgrid")

        ax = sns.relplot(data=self.interp_df[probe], kind='line')
        fig = plt.gcf()
        fig.set_size_inches(14, 8.66)
        xticks = ax.get_xticks()
        ax.set_xticklabels([pd.to_datetime(tm, unit='s').tz_localize('UTC').tz_convert('US/Central').strftime('%Y-%m-%d\n %H:%M:%S %Z')
                              for tm in xticks], rotation=0)
        return fig