#!/usr/bin/env python
# coding: utf-8

# # Analysis Helper

# In[ ]:


import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz

import gm2
import trfp
import helper_function_candidates as helper


# # Convert Tier-1 ROOT files to Pandas data frames

# This function takes a list of run numbers and uses S. Corrodi's gm2 module to read the ROOT files.
# It then performs the event-wise DQC cuts and time-grid interpolation.
# It outputs a standard Pandas data frame.

# ## Short helper functions
# `_nan_helper`
# 
# `_choose_J`

# In[ ]:


def _nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, index= nan_helper(y)
        >>> y[nans]= np.interp(index(nans), index(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def _choose_J(st):
    if trfp.STATION_PROBE_NUMBER[st] == 4:
        if st == 41: J = trfp.J_4_PROBE_ST41
        elif (st == 37) | (st == 39): J = trfp.J_4_PROBE_ST37_ST39
        else: J = trfp.J_4_PROBE
    elif trfp.STATION_PROBE_NUMBER[st] == 6:
        if st < 6: J = trfp.J_6_PROBE_OFFSET
        else: J = trfp.J_6_PROBE
    else:
        raise Exception('Invalid number of station probes.')
    return J

def _choose_theta(st):
    if trfp.STATION_PROBE_NUMBER[st] == 4:
        if st == 41:
            theta_fp = trfp.THETA_FP_4_ST41
        elif (st == 37) | (st == 39):
            theta_fp = trfp.THETA_FP_4_ST37_ST39
        else:
            theta_fp = trfp.THETA_FP_4
    elif trfp.STATION_PROBE_NUMBER[st] == 6:
        theta_fp = trfp.THETA_FP_6
    else:
        raise Exception('Invalid number of station probes.')
    return theta_fp
    


# ## Root to interpolated data frame
# `root_to_pandas`

# In[ ]:


def root_to_pandas(run_range, prefix=None, tr_run=False):
    
    if len(run_range) == 0: raise Exception('No runs specified.')
    if tr_run and len(run_range) > 1: raise Exception('Only one trolley run can be processed at a time.')
    if tr_run:
        tr_run = gm2.Trolley(run_range, prefix=prefix)
        if tr_run.getEntries() == 0: raise Exception('No trolley events.')
    fp_run = gm2.FixedProbe(run_range, prefix=prefix)
    
    if tr_run:
        tr_time, tr_phi, tr_freq = tr_run.getBasics()
        # drop first 5 events, which are 0, and last event, which can some times be 0
        tr_time = tr_time[5:-1,:]/1.0e9  # convert nsec to sec
        tr_phi = tr_phi[5:-1,:]
        tr_freq = tr_freq[5:-1,:]
        for tr in range(17):
            tr_freq[:, tr] += trfp.PLUNGING_PROBE_CALIBRATIONS[tr]

    fp_time, fp_freq, fp_qual = fp_run.getBasics()
    # drop first event, which is always 0 from import
    fp_time = fp_time[1:,:]/1.0e9  # convert nsec to sec
    fp_freq = fp_freq[1:,:]
    fp_qual = fp_qual[1:,:]

    ################################################################################
    ### Apply fixed probe event DQC cuts ###########################################
    ################################################################################
    
    # remove the 8th bit and 16th bit flags (Ran's flags?)
    fp_qual[fp_qual >= 2**16] -= 2**16
    fp_qual[fp_qual >= 2**8] -= 2**8

    fp_freq_dqc = fp_freq.copy()
    fp_freq_dqc[fp_qual > 0] = np.nan

    for fp in range(378):
        nans, index = _nan_helper(fp_freq_dqc[:,fp])
        fp_freq_dqc[:,fp][nans] = np.interp(index(nans), index(~nans), fp_freq_dqc[:,fp][~nans])

    ################################################################################
    ################################################################################
    ################################################################################


    # first, need to make array of raw times, grid times, and grid time edges
    # then find the index of each grid time point in that sorted array

    ################################################################################
    ### Trolley Runs ###############################################################
    ################################################################################
    if tr_run:
        grid_times = np.arange(np.ceil(np.max(tr_time[0,:])),
                               np.floor(np.min(tr_time[-1,:]))+1,
                               1.0)
        edge_times = np.arange(grid_times[0]-0.5, grid_times[-1]+1.5, 1.0)
        grid_phi = np.empty([grid_times.size, 1])
        grid_tr_freqs = np.empty([grid_times.size, 17])
        grid_freqs = np.empty([grid_times.size, 378])
   

        ################################################################################
        ### Interpolate trolley position ###############################################
        ################################################################################

        phi_time = np.mean(tr_time, axis=1)
        phi_phi = np.mean(tr_phi, axis=1)
        
        all_times = np.append(edge_times,phi_time)
        sort_index = np.argsort(all_times)
        unsort_index = np.argsort(sort_index)
        edge_index = unsort_index[0:len(grid_times)+1]
        
        edge_phi = np.interp(edge_times, phi_time, phi_phi)
        all_phi = np.append(edge_phi, phi_phi)
        sort_phi = all_phi[sort_index]
        sort_times = all_times[sort_index]
        
        integrated_phi = cumtrapz(sort_phi, x=sort_times, initial=0)
        grid_phi[:,0] = np.diff(integrated_phi[edge_index])

        ################################################################################
        ################################################################################
        ################################################################################
        
        for tr in range(17):
            all_times = np.append(edge_times, tr_time[:,tr])

            sort_index = np.argsort(all_times)
            unsort_index = np.argsort(sort_index)
            edge_index = unsort_index[0:len(grid_times)+1]

            # interpolate the freqs on the time grid edges
            edge_freqs = np.interp(edge_times, tr_time[:,tr], tr_freq[:,tr])

            # sort all times and all freqs by method calculated above
            all_freqs = np.append(edge_freqs, tr_freq[:,tr])
            sort_freqs = all_freqs[sort_index]
            sort_times = all_times[sort_index]

            # run cumtrapz on sorted times and freqs
            integrated_freqs = cumtrapz(sort_freqs, x=sort_times, initial=0)

            # take differences of the integrated freqs at grid edge points
            # to find the integral between the two points

            grid_tr_freqs[:,tr] = np.diff(integrated_freqs[edge_index])
        

        for fp in range(378):
            all_times = np.append(edge_times,fp_time[:,fp])

            sort_index = np.argsort(all_times)
            unsort_index = np.argsort(sort_index)
            edge_index = unsort_index[0:len(grid_times)+1]

            # interpolate the freqs on the time grid edges
            edge_freqs = np.interp(edge_times, fp_time[:,fp], fp_freq_dqc[:,fp])

            # sort all times and all freqs by method calculated above
            all_freqs = np.append(edge_freqs, fp_freq_dqc[:,fp])
            sort_freqs = all_freqs[sort_index]
            sort_times = all_times[sort_index]

            # run cumtrapz on sorted times and freqs
            integrated_freqs = cumtrapz(sort_freqs, x=sort_times, initial=0)

            # take differences of the integrated freqs at grid edge points
            # to find the integral between the two points

            grid_freqs[:,fp] = np.diff(integrated_freqs[edge_index])
        
        cols = ["tr_phi"] + ["tr" + str(i) for i in np.arange(17)] + ["fp" + str(i) for i in np.arange(378)]
        data = np.append(grid_phi, np.append(grid_tr_freqs, grid_freqs, axis=1), axis=1)
        
    ################################################################################
    ### Fixed Probe Runs ###########################################################
    ################################################################################
    else:
        grid_times = np.arange(np.ceil(np.max(fp_time[0,:])),
                               np.floor(np.min(fp_time[-1,:]))+1,
                               1.0)
        edge_times = np.arange(grid_times[0]-0.5, grid_times[-1]+1.5, 1.0)
        grid_freqs = np.empty([grid_times.size, 378])

        for fp in range(378):
            all_times = np.append(edge_times,fp_time[:,fp])

            sort_index = np.argsort(all_times)
            unsort_index = np.argsort(sort_index)
            edge_index = unsort_index[0:len(grid_times)+1]

            # interpolate the freqs on the time grid edges
            edge_freqs = np.interp(edge_times, fp_time[:,fp], fp_freq_dqc[:,fp])

            # sort all times and all freqs by method calculated above
            all_freqs = np.append(edge_freqs, fp_freq_dqc[:,fp])
            sort_freqs = all_freqs[sort_index]
            sort_times = all_times[sort_index]

            # run cumtrapz on sorted times and freqs
            integrated_freqs = cumtrapz(sort_freqs, x=sort_times, initial=0)

            # take differences of the integrated freqs at grid edge points
            # to find the integral between the two points

            grid_freqs[:,fp] = np.diff(integrated_freqs[edge_index])
            
        data = grid_freqs
        cols = ["fp" + str(i) for i in np.arange(378)]

    return pd.DataFrame(data, index=grid_times, columns=cols)


# ## Interpolated data frame to moment data frame
# `calc_moment_df`

# In[ ]:


def calc_moment_df(interp_df):
    
    tr_run = 'tr_phi' in interp_df.columns
    
    moment_df = pd.DataFrame(index=interp_df.index)
    
    # calculate trolley moments if needed
    if tr_run:
        moment_df['tr_phi'] = interp_df['tr_phi'].copy()
        print 'Calculating trolley moments.',
        theta_tr = trfp.THETA_TR
        for m in np.arange(17):
            tr_probes = ['tr'+str(probe) for probe in np.arange(17)]
            moment_df['tr,m'+str(m+1)] = interp_df[tr_probes].dot(theta_tr[m])

    # create the 72*6 fixed probe moments
    for station in np.arange(72):
        print '\rCalculating station ' + str(station) + ' moments.',
        fp_st = ['fp'+str(fp) for fp in trfp.STATION_PROBE_ID[station]]

        # choose proper theta matrix
        theta_fp = _choose_theta(station)

        # step through m values
        for m in np.arange(len(trfp.STATION_PROBE_ID[station])):
            stm = 'st'+str(station)+',m'+str(m+1)
            moment_df[stm] = interp_df[fp_st].dot(theta_fp[m])
        if len(trfp.STATION_PROBE_ID[station]) == 4:
            moment_df['st'+str(station)+',m5'] = np.nan
            moment_df['st'+str(station)+',m6'] = np.nan

    # Interpolate the 6-probe m5 and m6 into the 4-probe m5 and m6        
    for st in range(72):
        if trfp.STATION_PROBE_NUMBER[st] == 4:
            wt = trfp.STATION_BARCODE_PHI[(st+1)%72] - trfp.STATION_BARCODE_PHI[(st-1)%72]
            w1 = trfp.STATION_BARCODE_PHI[(st+1)%72] - trfp.STATION_BARCODE_PHI[st]
            w2 = trfp.STATION_BARCODE_PHI[st] - trfp.STATION_BARCODE_PHI[(st-1)%72]

            for m in [5, 6]:
                stm = 'st'+str(st)+',m'+str(m)
                if not np.isnan(moment_df[stm].iloc[0]):
                    print stm + ' is not nan.'
                    continue
                stm1 = 'st'+str((st-1)%72)+',m'+str(m)
                stm2 = 'st'+str((st+1)%72)+',m'+str(m)

                moment_df[stm] = (w1*moment_df[stm1] + w2*moment_df[stm2])/wt

    print '\rFinished calculating all moments for ' + str(moment_df.shape[0]) + ' events.'
    
    return moment_df


# # Virtual trolley measurement calculation

# `vtm_calc`
# 
# $m_{vtr}(t) = m_{tr}(0) + J\left[m_{fp}(t) - m_{fp}(0)\right]$
# 
# $m_{vtr}(t) = J m_{fp}(t) + \left[m_{tr}(0) - J m_{fp}(0)\right]$

# In[ ]:


# generate virutal trolley measurements

def vtm_calc(fp_moment_df,
             baseline_time_1, baseline_time_2,
             tr_baseline_1, tr_baseline_2,
             fp_baseline_1, fp_baseline_2):
    
    vtm_df = fp_moment_df.copy()
    
    # apply Jacobian to each station's m_fp
    
    for st in range(72):
        
        # choose Jacobian
        J = _choose_J(st)
                
        # apply Jacobian to each station's m_fp
        stms = ['st'+str(st)+',m'+str(m+1) for m in range(5)]
        vtm_df[stms] = vtm_df[stms].dot(np.transpose(J))
    
        # apply sync correction (m_tr(0) - J*m_fp(0)), interpolated
        
        sync_correction_1 = tr_baseline_1[st,0:5] - np.matmul(J, fp_baseline_1[st,0:5])
        sync_correction_2 = tr_baseline_2[st,0:5] - np.matmul(J, fp_baseline_2[st,0:5])
        
        sync_correction = np.empty((len(vtm_df.index.values),5))
        for m in range(5):
            sync_correction[:,m] = np.interp(vtm_df.index.values,
                                             [baseline_time_1[st,m], baseline_time_2[st,m]],
                                             [sync_correction_1[m], sync_correction_2[m]])
        vtm_df[stms] += sync_correction
        
    return vtm_df


# # Trolley footprint removal

# Remove the trolley footprint from the fixed probe meaurements by vetoing a window when the trolley is near the fixed probe station and replacing the vetoed section with a drift-corrected average from the rest of the ring.

# In[ ]:





# # Trolley run station averaging

# Functions that average trolley moments over azimuth and corrected fixed probe moments over time. Two seprate averaging modes: boxcar and triangular.
# 
# NOTE: more work is needed to determine how to applying triangular averaging to fixed probe time average.

# In[ ]:





# # Sync offset calculation

# In[ ]:


def sync_offset_calc(tr_corr_df_1, tr_corr_df_2):
    tr_baseline_1, fp_baseline_1, baseline_time_1, _, _ = helper.trolley_run_station_average(tr_corr_df_1)
    tr_baseline_2, fp_baseline_2, baseline_time_2, _, _ = helper.trolley_run_station_average(tr_corr_df_2)
    delta_time = baseline_time_2 - baseline_time_1
    
    delta_tr_baseline = tr_baseline_2 - tr_baseline_1
    delta_fp_baseline = fp_baseline_2 - fp_baseline_1
    
    sync_offsets = np.empty((72,5))
    for st in range(72):
        J = _choose_J(st)
#         sync_offsets[st,:] = delta_tr_baseline[st,:-1] - delta_fp_baseline[st,:-1]
        sync_offsets[st,:] = delta_tr_baseline[st,:-1] -  np.matmul(J, delta_fp_baseline[st,:-1])
    
    return sync_offsets, delta_time
