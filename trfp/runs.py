"""Contains classes for filling trolley and fixed probe runs."""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import scipy

import matplotlib.pyplot as plt
import seaborn as sns

import gm2
import trfp

from datetime import datetime
import pytz

# Super class for both trolley and fixed probe runs
class Run(object):
    
    """Class for building DataFrames of Trolley runs."""

    def __init__(self, run):
        self.run = run

        # determine whether the input run is a trolley run or fixed probe run
        tr_run_temp = gm2.Trolley([run])
        if tr_run_temp.getEntries() != 0:
            self.trolley = True
            self.tr_run = tr_run_temp
            print 'Trolley run.'
        else:
            self.trolley = False
            print 'Fixed probe run.'
        
        self.fp_run = gm2.FixedProbe([run])

        if self.trolley:
            times, tr_freq_interp, tr_phi_interp, fp_freq_interp = self.__time_interpolation_tr()
            # apply plunging probe calibrations to tr_freq_interp
            for ii in np.arange(17):
                tr_freq_interp[:,ii] = tr_freq_interp[:,ii] + trfp.PLUNGING_PROBE_CALIBRATIONS[ii]
            cols = ["tr_phi"] + ["tr" + str(i) for i in np.arange(17)] + ["fp" + str(i) for i in np.arange(378)]
            data = np.append(np.append(tr_phi_interp, tr_freq_interp, axis=1),
                             fp_freq_interp, axis=1)
            self.interp_df = pd.DataFrame(data, index=times, columns=cols)
        else:
            times, fp_freq_interp = self.__time_interpolation_fp()
            cols = ["fp" + str(i) for i in np.arange(378)]
            data = fp_freq_interp
            self.interp_df = pd.DataFrame(data, index=times, columns=cols)
            
#         # add uncertainty column to interp_df for each column
#         for column in list(self.interp_df):
#             self.interp_df[column+'_err'] = np.nan
        
        self.moment_df = self.__moment_dataframe()

    def __time_interpolation_tr(self):
        tr_time, tr_phi, tr_freq = self.tr_run.getBasics(mode_phi=2)
        fp_time, fp_freq, fp_qual = self.fp_run.getBasics()
        tr_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds
        fp_time /= 1.0e9
        
        fp_qual[fp_qual == 2**16] = 0
        fp_qual[fp_qual == 2**8] = 0

        fp_freq[fp_qual > 0] = np.nan  # veto bad fids

        tr_indices = np.arange(len(tr_freq)) >= 10  # drop first 10 trolley events
        fp_indices = np.arange(len(fp_freq)) >= 3  # drop first 3 fixed probe events

        times = np.arange(np.ceil(np.max([tr_time[tr_indices][0, 16],
                                          fp_time[fp_indices][0, 377]])),
                          np.floor(np.min([tr_time[tr_indices][-1, 0],
                                           fp_time[fp_indices][-1, 0]])) + 1,
                          1)
        steps = 100
        rate = 1./steps
        offset = steps/2

        integration_times = np.arange(times[0], times[-1], rate)

        tr_freq_interp = np.zeros([times.size, 17])
        tr_phi_interp = np.zeros([times.size, 1])
        fp_freq_interp = np.zeros([times.size, 378])
        
        # integrate trolley phi measurements
        # unwrap trolley measurement
        print 'Interpolating trolley position.'
        index = np.arange(tr_phi[1:,0].size)[np.abs(tr_phi[1:,0]-tr_phi[0:-1,0]) > 300]
        if index.size == 1:
            index = int(index)
            tr_phi_out = tr_phi[:,0].copy()
            tr_phi_out[index+1:] = tr_phi_out[index+1:] + np.sign(tr_phi_out[index]-tr_phi_out[index+1])*360
        elif index.size == 0:
            tr_phi_out = tr_phi[:,0].copy()
        else:
            raise TrolleyPositionError('More than one wrap around.')
        
        phi_interp = scipy.interpolate.interp1d(tr_time[:,0], tr_phi_out, kind='slinear')
        integration_points = phi_interp(integration_times)

        cumulative_integration = scipy.integrate.cumtrapz(integration_points, integration_times,
                                                         initial=0) + integration_points[0]

        indices = [steps*i + offset for i in range(integration_times.size/steps)]

        output = np.diff(cumulative_integration[indices])
        output = np.append(np.array([integration_points[0]]), output)
        output = np.append(output, np.array([integration_points[-1]]))

        tr_phi_interp[:,0] = output
        tr_phi_interp = tr_phi_interp%360
        
        # integrate each trolley probe
        print 'Interpolating trolley frequencies.'
        for probe in range(17):
            probe_interp = scipy.interpolate.interp1d(tr_time[:,probe], tr_freq[:,probe], kind='cubic')
            integration_points = probe_interp(integration_times)

            cumulative_integration = scipy.integrate.cumtrapz(integration_points, integration_times,
                                                             initial=0) + integration_points[0]

            indices = [steps*i + offset for i in range(integration_times.size/steps)]

            output = np.diff(cumulative_integration[indices])
            output = np.append(np.array([integration_points[0]]), output)
            output = np.append(output, np.array([integration_points[-1]]))

            tr_freq_interp[:,probe] = output
        
        # integrate each fixed probe
        print 'Interpolating fixed probe frequencies.'
        for probe in range(378):
            not_nan = ~np.isnan(fp_freq[:,probe])
            integration_points = np.interp(integration_times, fp_time[:,probe][not_nan], fp_freq[:,probe][not_nan])

            cumulative_integration = scipy.integrate.cumtrapz(integration_points, integration_times,
                                                             initial=0) + integration_points[0]

            indices = [steps*i + offset for i in range(integration_times.size/steps)]

            output = np.diff(cumulative_integration[indices])
            output = np.append(np.array([integration_points[0]]), output)
            output = np.append(output, np.array([integration_points[-1]]))

            fp_freq_interp[:,probe] = output

        return times, tr_freq_interp, tr_phi_interp, fp_freq_interp
    
    def __time_interpolation_fp(self):
        fp_time, fp_freq, fp_qual = self.fp_run.getBasics()
        fp_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds
        
        fp_qual[fp_qual == 2**16] = 0
        fp_qual[fp_qual == 2**8] = 0
        
        fp_freq[fp_qual > 0] = np.nan  # veto bad fids

        fp_indices = np.arange(len(fp_freq)) >= 1  # drop first 3 fixed probe events

        times = np.arange(np.ceil(fp_time[fp_indices][0, 377]),
                          np.floor(fp_time[fp_indices][-1, 0]) + 1,
                          1)  # NOTE THIS IS 1 SECOND NOW, DOESN'T MATTER BECAUSE OF INTEGRATION METHOD
        steps = 100
        rate = 1./steps
        offset = steps/2

        integration_times = np.arange(times[0], times[-1], rate)

        fp_freq_interp = np.empty([times.size, 378])
        print 'Interpolating fixed probe frequencies.'
        
        for probe in range(378):
            not_nan = ~np.isnan(fp_freq[:,probe])
            integration_points = np.interp(integration_times, fp_time[:,probe][not_nan], fp_freq[:,probe][not_nan])

            cumulative_integration = scipy.integrate.cumtrapz(integration_points, integration_times,
                                                             initial=0) + integration_points[0]

            indices = [steps*i + offset for i in range(integration_times.size/steps)]

            output = np.diff(cumulative_integration[indices])
            output = np.append(np.array([integration_points[0]]), output)
            output = np.append(output, np.array([integration_points[-1]]))

            fp_freq_interp[:,probe] = output
            
        return times, fp_freq_interp


    def __moment_dataframe(self):
        """Builds dataframe of field moments."""

        moment_df = pd.DataFrame(index=self.interp_df.index)
        
        if self.trolley:
            # copy the interpolation data for the phi
            tr_phi = self.interp_df['tr_phi'].copy()
            moment_df['tr_phi'] = self.interp_df['tr_phi'].copy()
#             moment_df = pd.DataFrame(tr_phi, index=tr_phi.index, columns=['tr_phi'])

            # create the 17 trolley moments at each point in time
            print 'Calculating trolley moments.',
            theta_tr = trfp.THETA_TR
            for m in np.arange(17):

                tr_probes = ['tr'+str(probe) for probe in np.arange(17)]
                moment_df['tr,m'+str(m+1)] = self.interp_df[tr_probes].dot(theta_tr[m])

        # create the 72*6 fixed probe moments
        for station in np.arange(72):
            print '\rCalculating station ' + str(station) + ' moments.',
            fp_st = ['fp'+str(fp) for fp in trfp.STATION_PROBE_ID[station]]
            
            # choose proper theta matrix
            if len(trfp.STATION_PROBE_ID[station]) == 4:
                if station == 41:
                    theta_fp = trfp.THETA_FP_4_ST41
                elif (station == 37) | (station == 39):
                    theta_fp = trfp.THETA_FP_4_ST37_ST39
                else:
                    theta_fp = trfp.THETA_FP_4
            else:
                theta_fp = trfp.THETA_FP_6
                
            # step through m values
            for m in np.arange(len(trfp.STATION_PROBE_ID[station])):
                stm = 'st'+str(station)+',m'+str(m+1)
                moment_df[stm] = self.interp_df[fp_st].dot(theta_fp[m])
            if len(trfp.STATION_PROBE_ID[station]) == 4:
                moment_df['st'+str(station)+',m5'] = np.nan
                moment_df['st'+str(station)+',m6'] = np.nan

        print '\rFinished calculating all moments for ' + str(moment_df.shape[0]) + ' events.'

        return moment_df

    def save_h5(self, file_name):
        """Save interpolation data frame and moment dataframe in hdf5 format."""
        interp_key = "run_" + str(self.run) + "_interp_df"
        moment_key = "run_" + str(self.run) + "_moment_df"
        
        print "Saving run " + str(self.run) + '.'
        
        self.interp_df.to_hdf(file_name, key=interp_key)
        self.moment_df.to_hdf(file_name, key=moment_key)
        
#     # The following are slower methods for time-grid interpolation
    
#     def __time_interpolation_tr(self):
#         tr_time, tr_phi, tr_freq = self.tr_run.getBasics(mode_phi=2)
#         _, fp_time, fp_freq = self.fp_run.getBasics()
#         tr_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds
#         fp_time /= 1.0e9

#         tr_indices = np.arange(len(tr_freq)) >= 10  # drop first 10 trolley events
#         fp_indices = np.arange(len(fp_freq)) >= 3  # drop first 3 fixed probe events

#         times = np.arange(np.ceil(np.max([tr_time[tr_indices][0, 16],
#                                           fp_time[fp_indices][0, 377]])),
#                           np.floor(np.min([fp_time[fp_indices][-1, 0],
#                                            fp_time[fp_indices][-1, 0]])) + 1,
#                           1)
#         time_edges = (times + np.roll(times,1))/2
#         time_edges[0] = times[0]
#         time_edges = np.append(time_edges, times[-1])
#         delta_ts = time_edges[1:] - time_edges[0:-1]

#         tr_freq_interp = np.zeros([times.size, 17])
#         tr_phi_interp = np.zeros([times.size, 1])
#         fp_freq_interp = np.zeros([times.size, 378])
        
#         # integrate trolley phi measurements
#         print 'Interpolating trolley phi.',
#         phi_interp = scipy.interpolate.interp1d(tr_time[:,0], tr_phi[:,0], kind='cubic')
#         for t in np.arange(times.shape[0]):
#             tr_phi_interp[t,0], _ = scipy.integrate.fixed_quad(phi_interp, time_edges[t], time_edges[t+1], n=50)
#             # n=50 for ppb level accuracy
#         tr_phi_interp[:,0] /= delta_ts
#         print '\n'
#         tr_phi_interp = tr_phi_interp%360
        
#         # integrate each trolley probe
#         for tr in np.arange(17):
#             print '\rInterpolating trolley probe '+str(tr)+ '.',
#             probe_interp = scipy.interpolate.interp1d(tr_time[:,tr], tr_freq[:,tr], kind='cubic')
#             for t in np.arange(times.shape[0]):
#                 tr_freq_interp[t, tr], _ = scipy.integrate.fixed_quad(probe_interp, time_edges[t], time_edges[t+1], n=50)
#                 # n=50 for ppb level accuracy
#             tr_freq_interp[:,tr] /= delta_ts
#         print '\n'
        
#         # integrate each fixed probe
#         for fp in np.arange(378):
#             print '\rInterpolating fixed probe '+str(fp)+ '.',
#             probe_interp = scipy.interpolate.interp1d(fp_time[:,fp], fp_freq[:,fp], kind='cubic')
#             for t in np.arange(times.shape[0]):
#                 fp_freq_interp[t, fp], _ = scipy.integrate.fixed_quad(probe_interp, time_edges[t], time_edges[t+1], n=50)
#                 # n=50 for ppb level accuracy
#             fp_freq_interp[:,fp] /= delta_ts
#         print '\n'

#         return times, tr_freq_interp, tr_phi_interp, fp_freq_interp
    
#     def __time_interpolation_fp(self):
#         _, fp_time, fp_freq = self.fp_run.getBasics()
#         fp_time /= 1.0e9  # timestamps come in nanoseconds, convert to seconds

#         fp_indices = np.mean(fp_freq, 1) > 0

#         times = np.arange(np.ceil(fp_time[fp_indices][0, 377]),
#                           np.floor(fp_time[fp_indices][-1, 0]) + 1,
#                           2)
#         time_edges = (times + np.roll(times,1))/2
#         time_edges[0] = times[0]
#         time_edges = np.append(time_edges, times[-1])
#         delta_ts = time_edges[1:] - time_edges[0:-1]

#         fp_freq_interp = np.zeros([times.size, 378])

#         # integrate each fixed probe
#         for fp in np.arange(378):
#             print '\rInterpolating fixed probe '+str(fp)+ '.',
#             probe_interp = scipy.interpolate.interp1d(fp_time[:,fp], fp_freq[:,fp], kind='cubic')
#             for t in np.arange(times.shape[0]):
#                 fp_freq_interp[t, fp], _ = scipy.integrate.fixed_quad(probe_interp, time_edges[t], time_edges[t+1], n=50)
#                 # n=50 for ppb level accuracy
#             fp_freq_interp[:,fp] /= delta_ts
#         print '\n'

#         return times, fp_freq_interp
    
