import gm2
import trfp
import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz
import psycopg2

###
### Short utility functions
###

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
        if st < 6:
            if st==5: J = trfp.J_ST_5
            else: J = trfp.J_6_PROBE_OFFSET
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

def _get_keys(filename):
    hdf = pd.HDFStore(filename, mode='r')
    keys = hdf.keys()
    hdf.close()
    output_keys = []
    for key in keys:
        if key != '/subrun_df': output_keys.append(key[1:])
    return output_keys

###
### Functions that convert ROOT data products in Python (pandas) data products
###

def root_to_pandas(run_range, prefix=None, tr_run=False, sanitize=False):
    
    if len(run_range) == 0: raise Exception('No runs specified.')
    if tr_run and len(run_range) > 1: raise Exception('Only one trolley run can be processed at a time.')
    if tr_run:
        tr_run = gm2.Trolley(run_range, prefix=prefix)
        if tr_run.getEntries() == 0: raise Exception('No trolley events.')
    fp_run = gm2.FixedProbe(run_range, prefix=prefix)
    
    if tr_run:
        # drop first 5 events, which are 0, and last event, which can some times be 0
        tr_time = tr_run.time[5:-1,:]/1.0e9  # convert nsec to sec
        tr_phi = tr_run.azi[5:-1,:]
        tr_freq = tr_run.freq[5:-1,:]
        for tr in range(17):
            tr_freq[:, tr] += trfp.PLUNGING_PROBE_CALIBRATIONS[tr]

#     fp_time, fp_freq, fp_qual = fp_run.getBasics()
    # drop first event, which is always 0 from import
    fp_time = fp_run.time[1:,:]/1.0e9  # convert nsec to sec
    fp_freq = fp_run.freq[1:,:]
    fp_qual = fp_run.qtag[1:,:]

    ################################################################################
    ### Apply fixed probe event DQC cuts ###########################################
    ################################################################################
    
    # remove the 8th bit and 16th bit flags (Ran's flags?)
    fp_qual[fp_qual >= 2**16] -= 2**16
    fp_qual[fp_qual >= 2**8] -= 2**8

    fp_freq_dqc = fp_freq.copy()
    fp_freq_dqc[fp_qual > 0] = np.nan
    
    ################################################################################
    ### Apply fixed probe event hell probe cuts ####################################
    ### Interpolate over nan values ################################################
    ################################################################################
    
    for fp in range(378):
        if sanitize:
            fp_freq_dqc[:,fp], _ = _drop_diff_single(fp_freq_dqc[:,fp], threshold=200)
            fp_freq_dqc[:,fp], _ = _drop_freq_bin(fp_freq_dqc[:,fp], bin_len = 1000, std = 3)
        
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

        phi_time = np.median(tr_time, axis=1)
        phi_phi = np.median(tr_phi, axis=1)
        
        pivot = np.argmax(np.abs(np.diff(phi_phi)))
        phi_phi[pivot+1:] -= 360
        
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
        
        grid_phi[grid_phi < 0] += 360

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

###
### Averaging functions
###

def trolley_run_station_average(corrected_df):
        
    station_edges = trfp.STATION_BARCODE_EDGES
    station_edges_6_probe = trfp.STATION_BARCODE_EDGES_6

    # tr_phi is not monotonic, so sort by tr_phi

    corrected_df = corrected_df.sort_values(by=['tr_phi'])

    measured_phi = corrected_df['tr_phi'].values
    measured_extent = (np.roll(measured_phi,-1)-np.roll(measured_phi,1))/2
    measured_extent[0] = measured_extent[0]+180
    measured_extent[-1] = measured_extent[-1]+180
    # print np.max(measured_extent)

    corrected_df['tr_extent'] = pd.Series(measured_extent, index=corrected_df.index)
    corrected_df = corrected_df.sort_index()

    # for a given station:
    # create a mask for when trolley is in [low edge, high edge)
    tr_baseline = np.full((72,9), np.nan)
    fp_baseline = np.full((72,5), np.nan)
    summed_azimuth = np.full(72, np.nan)
    summed_pts = np.full(72, np.nan)
    baseline_time = np.full(72, np.nan)

    st6 = 0
    for st in range(72):
        
        num_probes = len(trfp.STATION_PROBE_ID[st])
        # first do m1-4 for all stations
        
        if station_edges[st+1] > station_edges[st]:
            mask = (corrected_df['tr_phi'] >= station_edges[st]) & (corrected_df['tr_phi'] < station_edges[st+1])
        else:  # case where we go over the 360 deg line
            mask = (corrected_df['tr_phi'] >= station_edges[st]) | (corrected_df['tr_phi'] < station_edges[st+1])

        out_df = corrected_df[mask].copy()
        summed_pts[st] = out_df.shape[0]
        summed_azimuth[st] = sum(out_df['tr_extent'].values)
        baseline_time[st] = sum(out_df.index.values)/summed_pts[st]

        # not using fp m6
        for m in range(5):
            st_id = 'tr,m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                tr_baseline[st, m] = sum(out_df['tr_extent'].values*out_df[st_id].values)/sum(out_df['tr_extent'].values)

            st_id = 'st'+str(st)+',m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                fp_baseline[st, m] = np.mean(out_df[st_id])
                
        for m in range(5,9):
            st_id = 'tr,m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                tr_baseline[st, m] = sum(out_df['tr_extent'].values*out_df[st_id].values)/sum(out_df['tr_extent'].values)

    return tr_baseline, fp_baseline, baseline_time, summed_azimuth, summed_pts

###
### Dataframe manipulation
###

def _apply_blinds_fp(input_df, blinds):
    output_df = input_df.copy()
    for m in range(6):
        stms = ['st'+str(st)+',m'+str(m+1) for st in range(72)]
        output_df[stms] = output_df[stms] + blinds[m]
    return output_df

def _apply_blinds_tr(input_df, blinds):
    output_df = input_df.copy()
    for m in range(6):
        stms = ['st'+str(st)+',m'+str(m+1) for st in range(72)]
        output_df[stms] = output_df[stms] + blinds[m]
        trms = ['tr,m'+str(m+1)]
        output_df[trms] = output_df[trms] + blinds[m]
    return output_df

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
        # note that this no longer calculates m6 for 6 probe stations
        ### CLEAN THIS PART UP
        for m in np.arange(len(trfp.STATION_PROBE_ID[station])):
            if not m == 5:
                stm = 'st'+str(station)+',m'+str(m+1)
                moment_df[stm] = interp_df[fp_st].dot(theta_fp[m])
        if len(trfp.STATION_PROBE_ID[station]) == 4:
            moment_df['st'+str(station)+',m5'] = np.nan

    # Interpolate the 6-probe m5 and m6 into the 4-probe m5 and m6        
    for st in range(72):
        if trfp.STATION_PROBE_NUMBER[st] == 4:
            wt = trfp.STATION_BARCODE_PHI[(st+1)%72] - trfp.STATION_BARCODE_PHI[(st-1)%72]
            w1 = trfp.STATION_BARCODE_PHI[(st+1)%72] - trfp.STATION_BARCODE_PHI[st]
            w2 = trfp.STATION_BARCODE_PHI[st] - trfp.STATION_BARCODE_PHI[(st-1)%72]

            # again, no m6 in fixed probes
            for m in [5]:
                stm = 'st'+str(st)+',m'+str(m)
                if not np.isnan(moment_df[stm].iloc[0]):
                    print stm + ' is not nan.'
                    continue
                stm1 = 'st'+str((st-1)%72)+',m'+str(m)
                stm2 = 'st'+str((st+1)%72)+',m'+str(m)

                moment_df[stm] = (w1*moment_df[stm1] + w2*moment_df[stm2])/wt

    print '\rFinished calculating all moments for ' + str(moment_df.shape[0]) + ' events.'
    
    return moment_df

def trolley_footprint_replacement(tr_moment_df, veto_extent=25):
    no_mask_df = tr_moment_df.copy()
    t0 = no_mask_df.index.values[0]
    no_mask_df.index -= t0
    index = no_mask_df.index.values

    for st in range(72):
        print '\rRemoving trolley image from station '+str(st)+'.',

        # veto when trolley is close to station

        veto_low = (trfp.STATION_BARCODE_PHI[st]-1.5-veto_extent/2)%360
        veto_high = (trfp.STATION_BARCODE_PHI[st]-1.5+veto_extent/2)%360
        if veto_low < veto_high:
            veto_mask = (no_mask_df['tr_phi']>veto_low) & (no_mask_df['tr_phi']<veto_high)
        else:  # this happens when wrapping around 360 deg
            veto_mask = (no_mask_df['tr_phi']>veto_low) | (no_mask_df['tr_phi']<veto_high)

        # no longer dealing with m6 in the fixed probes
        for m in range(5):

            stm = 'st'+str(st)+',m'+str(m+1)

            # calculate local drift

            times = no_mask_df.index.values[~veto_mask]
            freqs = no_mask_df[stm][~veto_mask]

            local_drift_fit = np.polyfit(times, freqs, 5)
            local_drift = np.polyval(local_drift_fit, no_mask_df.index.values)

            # need to average other side of ring
            all_good_stations = np.arange(6,72)  # not using the inflector stations
            no_ground_loop_stations = np.array(range(6,16)+range(64,72))  # vaid for 25 deg veto

            # next need to average all good stations that are not within 3 of current station
            if st not in range(16, 23):  # note that these ranged were chosen for 25 deg veto
                averaging_stations = np.delete(all_good_stations,
                                               np.argwhere((np.abs((all_good_stations - st)%72)<=3)
                                                          | (np.abs((all_good_stations - st)%72)>=69))
                                              )
            else:
                averaging_stations = np.delete(no_ground_loop_stations,
                                               np.argwhere((np.abs((no_ground_loop_stations - st)%72)<=3)
                                                          | (np.abs((no_ground_loop_stations - st)%72)>=69))
                                              )
            avg_stms = ['st'+str(avg_st)+',m'+str(m+1) for avg_st in averaging_stations]
            replacement = no_mask_df[avg_stms].mean(axis=1)

            # calculate global drift
            global_drift_fit = np.polyfit(index[veto_mask], replacement[veto_mask], 1)
            global_drift = np.polyval(global_drift_fit, index[veto_mask])

            # subtract global drift from replacement
            replacement = replacement[veto_mask] - global_drift

            # add local drift
            replacement = replacement + local_drift[veto_mask]

            no_mask_df[stm][veto_mask] = replacement

    no_mask_df.index += t0
    return no_mask_df


def vtm_calc(fp_moment_df,
             baseline_time_1, baseline_time_2,
             tr_baseline_1, tr_baseline_2,
             fp_baseline_1, fp_baseline_2):
    
    vtm_df = fp_moment_df.copy()
    
    # apply Jacobian to each station's m_fp
    
    # this is for moments 1 through 5 (fp trackable moments)
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
                                             [baseline_time_1[st], baseline_time_2[st]],
                                             [sync_correction_1[m], sync_correction_2[m]])
        vtm_df[stms] += sync_correction
    
    # linear interpolation of the trolley baselines for fp-untrackable moments
    for st in range(72):
        for m in range(5,9):
            stm = 'st'+str(st)+',m'+str(m+1)
            high_m_interp = np.interp(vtm_df.index.values, [baseline_time_1[st], baseline_time_2[st]],
                                      [tr_baseline_1[st, m], tr_baseline_2[st, m]])
            vtm_df[stm] = high_m_interp
    
    vtm_df = vtm_df.reindex(sorted(vtm_df.columns), axis=1)
        
    return vtm_df

def bloch_style_moments(corrected_dfs, keys):
    print '\nImplementing Bloch-style treatment of stations 1, 3, and 5.'
    station_phi = trfp.STATION_BARCODE_PHI
    weight10 = (station_phi[2]-station_phi[1])/((station_phi[2]-station_phi[0]))
    weight12 = (station_phi[1]-station_phi[0])/((station_phi[2]-station_phi[0]))
    weight32 = (station_phi[4]-station_phi[3])/((station_phi[4]-(station_phi[2]-360)))
    weight34 = (station_phi[3]-(station_phi[2]-360))/((station_phi[4]-(station_phi[2]-360)))
    weight54 = (station_phi[6]-station_phi[5])/((station_phi[6]-station_phi[4]))
    weight56 = (station_phi[5]-station_phi[4])/((station_phi[6]-station_phi[4]))
    
    bloch_style_dfs = {}

    for key in keys:
        print key
        bloch_style_dfs[key] = corrected_dfs[key].copy()
        
        # not using fp m6
        for m in range(1,6):
            print 'm' +str(m)+'\r',
            bloch_style_dfs[key]['st1,m'+str(m)] = weight10*corrected_dfs[key]['st0,m'+str(m)]+ weight12*corrected_dfs[key]['st2,m'+str(m)]
            bloch_style_dfs[key]['st3,m'+str(m)] = weight32*corrected_dfs[key]['st2,m'+str(m)] + weight34*corrected_dfs[key]['st4,m'+str(m)]
            bloch_style_dfs[key]['st5,m'+str(m)] = weight54*corrected_dfs[key]['st4,m'+str(m)] + weight56*corrected_dfs[key]['st6,m'+str(m)]
#             bloch_style_dfs[key]['st54,m'+str(m)] = weight5453*corrected_dfs[key]['st53,m'+str(m)] + weight5455*corrected_dfs[key]['st55,m'+str(m)]
    
    return bloch_style_dfs

###
### Data quality functions
###

def _drop_diff_single(input_array, threshold=200):
    #copy input
    start_array = input_array.copy()
    #sanitize input array by replacing nan values with median values
    start_array[np.isnan(start_array)] = np.nanmedian(start_array)
    diff = np.diff(start_array)
    drop_list = []
    
    for i in range(len(diff)-1):
        if np.abs(diff[i]) > threshold:
            drop_list.append(i+1)
    
    output_array = input_array.copy()
    output_array[drop_list] = np.nan
    
    return output_array, len(drop_list)

def _drop_freq_bin(input_array, bin_len = 1000, std = 3):
    start_array = input_array.copy()
    num_bin = start_array.size/bin_len

    for i in range(num_bin):

        center = np.nanmedian(start_array[i*bin_len:(i+1)*bin_len])
        width = np.nanstd(start_array[i*bin_len:(i+1)*bin_len])

        drop = np.abs(start_array[i*bin_len:(i+1)*bin_len] - center) > std*width
        start_array[i*bin_len:(i+1)*bin_len][drop] = np.nan

    center = np.nanmedian(start_array[num_bin*bin_len:])
    width = np.nanstd(start_array[num_bin*bin_len:])

    drop = np.abs(start_array[num_bin*bin_len:] - center) > std*width
    start_array[num_bin*bin_len:][drop] = np.nan
    
    drops = np.sum(np.isnan(start_array)) - np.sum(np.isnan(input_array))
    
    return start_array, drops

###
### Uncertainty functions
###

def sync_offset_calc(tr_corr_df_1, tr_corr_df_2):
    tr_baseline_1, fp_baseline_1, baseline_time_1, _, _ = trolley_run_station_average(tr_corr_df_1)
    tr_baseline_2, fp_baseline_2, baseline_time_2, _, _ = trolley_run_station_average(tr_corr_df_2)
    delta_time = baseline_time_2 - baseline_time_1
    
    delta_tr_baseline = tr_baseline_2 - tr_baseline_1
    delta_fp_baseline = fp_baseline_2 - fp_baseline_1
    
    sync_offsets = np.empty((72,9))
    for st in range(72):
        J = _choose_J(st)
        sync_offsets[st,0:5] = delta_tr_baseline[st,0:5] -  np.matmul(J, delta_fp_baseline[st,0:5])
        sync_offsets[st,5:9] = delta_tr_baseline[st,5:9]
    
    return sync_offsets, delta_time

###
### Import/export functions
###

def read_dfs(filename):
    keys = _get_keys(filename)
    interp_dfs = {}
    for key in keys:
        print '\nReading ' + key
        interp_dfs[key] = pd.read_hdf(filename, key=key)
    subrun_df = pd.read_hdf(filename, key='subrun_df')
    return interp_dfs, keys, subrun_df

def save_df(dict_run, save_dir, save_name, prefix='FieldPlainRootOutput_', sanitize=True):
    for frame in dict_run:
        print frame
        if frame != 'subrun_df':
            if frame[:2] == 'tr':
                new_df = root_to_pandas(dict_run[frame], prefix=prefix, tr_run=True, sanitize=sanitize)
            elif frame[:2] == 'fp':
                new_df = root_to_pandas(dict_run[frame], prefix=prefix, tr_run=False, sanitize=sanitize)
            else:
                raise NameError('Unexpected key name.')
        else:
            new_df = get_subrun_df(dict_run[frame][0], dict_run[frame][1])
        
        filename = save_dir+save_name
        new_df.to_hdf(filename, key=frame)
        
        print 'Saved to: '+filename

###
### High level wrappers
###

def interp_to_moment(interp_dfs, keys):
    moment_dfs = {}
    for key in keys:
        print '\nCalculating moments for ' + key
        moment_dfs[key] = calc_moment_df(interp_dfs[key])
    return moment_dfs

def blind_moments(moment_dfs, keys, blinds):
    blinded_moment_dfs = {}
    for key in keys:
        print '\nBlinding ' + key
        if key[:2] == 'tr': blinded_moment_dfs[key] = _apply_blinds_tr(moment_dfs[key], blinds)
        elif key[:2] == 'fp': blinded_moment_dfs[key] = _apply_blinds_fp(moment_dfs[key], blinds)
        else: raise NameError('Unexpected key name.')
    return blinded_moment_dfs
            
def moment_to_corrected(moment_dfs, keys):
    corrected_dfs = {}
    for key in keys:
        if key[:2] == 'tr': 
            print "\nRemoving trolley footprints for " + key
            corrected_dfs[key] = trolley_footprint_replacement(moment_dfs[key])
        elif key[:2] == 'fp': corrected_dfs[key] = moment_dfs[key].copy()
    print '\n'
    return corrected_dfs

def station_average(corrected_dfs, keys):
    '''
    Ultimately, the baselines should be stored in a dataframe with the old array names as rows.
    This will require rewriting `helper.vtm_calc` to use dfs instead of arrays.
    This will help keep everything in the clean dictionary of dfs framework.
    
    Could also update `helper_old.trolley_run_station_average`.
    '''
    
    print '\nCalculating trolley run baselines.'
    
    tr_baselines = {}
    fp_baselines = {}
    baseline_times = {}
    summed_azimuths = {}
    summed_pts = {}
    
#     ## future-proofing
#     station_avg_dfs = {}
    
    for key in keys:
        if key[:2] == 'tr':
            tr_baselines[key], fp_baselines[key], baseline_times[key], summed_azimuths[key], summed_pts[key] = trolley_run_station_average(corrected_dfs[key])

    baselines = {'tr':tr_baselines, 'fp':fp_baselines, 'time':baseline_times, 'azi':summed_azimuths, 'pts':summed_pts}
            
    return baselines

def calculate_vtms(corrected_dfs, keys, baselines, pair_dict):
    '''
    Going to need some new dictionaries that define run pairs (and single-sided runs).
    NOTE THIS NEEDS UPDATING FOR SINGLE SIDED RUNS.
    Future update: Make this work on both fixed probe and trolley runs (might need to update `helper.vtm_calc`).
    '''
    
    print '\nCalculating VTMs'
    
    vtm_dfs = {}
    
    for key in keys:
        if key[:2] == 'fp':
            vtm_dfs[key] = vtm_calc(corrected_dfs[key],
                                    baselines['time'][pair_dict[key][0]], baselines['time'][pair_dict[key][1]],
                                    baselines['tr'][pair_dict[key][0]], baselines['tr'][pair_dict[key][1]],
                                    baselines['fp'][pair_dict[key][0]], baselines['fp'][pair_dict[key][1]])
    
    return vtm_dfs

###
### Non-field functions
###

def get_subrun_df(start_run, end_run):
    
    dsn  = "dbname=gm2_online_prod user=gm2_reader host=g2db-priv port=5433"
    conn = psycopg2.connect(dsn)
    curr = conn.cursor()
    
    # get times by subrun
    subrun_time_columns = ['run', 'subrun', 'start_time', 'end_time', 'start_gps', 'end_gps']
    sql = "select "+", ".join(subrun_time_columns)+" from gm2dq.subrun_time where run >= %i and run <= %i order by run, subrun" % (start_run, end_run)
    curr.execute(sql)
    conn.commit()
    subrun_time = curr.fetchall()
    subrun_time_df = pd.DataFrame(subrun_time, columns=subrun_time_columns)
    
    # get ctags by subruns
    ctagswithdqc_columns = ['run', 'subrun', 'ctags', 't0val', 'fills']
    sql = "select "+", ".join(ctagswithdqc_columns)+" from gm2dq.ctagswithdqc where run >= %i and run <= %i order by run, subrun" % (start_run, end_run)
    curr.execute(sql)
    conn.commit()
    ctagswithdqc = curr.fetchall()
    ctagswithdqc_df = pd.DataFrame(ctagswithdqc, columns=ctagswithdqc_columns)
    
    # get subrun status database
    subrun_status_columns = ['run', 'subrun', 'quad_condition', 'kicker_condition', 'quad_ok',
                             'ctags_ok', 'losses_ok', 'fillcuts_ok', 'field_ok', 'trolley_period', 'field_period',
                             'ctags_loose_ok', 'quad_loose_ok', 'ctags_repeat_ok', 'losses_repeat_ok', 'fillcuts_repeat_ok']
    sql = "select "+", ".join(subrun_status_columns)+" from gm2dq.subrun_status where run >= %i and run <= %i order by run, subrun" % (start_run, end_run)
    curr.execute(sql)
    conn.commit()
    subrun_status = curr.fetchall()
    subrun_status_df = pd.DataFrame(subrun_status, columns=subrun_status_columns)
    
    # merge times, ctags, status into one subrun dataframe
    subrun_df = pd.merge(subrun_time_df, ctagswithdqc_df)
    subrun_df = pd.merge(subrun_df, subrun_status_df)
    subrun_df['ok'] = subrun_df['quad_ok'] & subrun_df['ctags_ok'] & subrun_df['losses_ok'] & subrun_df['fillcuts_ok'] & subrun_df['field_ok']
    subrun_df['start_time'] = subrun_df['start_time'].astype(np.int64)/1e9 + 5*60*60
    subrun_df['end_time'] = subrun_df['end_time'].astype(np.int64)/1e9 + 5*60*60
    subrun_df['start_gps'] = subrun_df['start_gps'].astype(np.int64)/1.0e9 + 5*60*60
    subrun_df['end_gps'] = subrun_df['end_gps'].astype(np.int64)/1.0e9 + 5*60*60
    
    print subrun_df.shape
    
    return subrun_df
