import numpy as np
import scipy

import pandas as pd

import trfp

# import functions

def import_fp_moment_runs(file_path, all_runs, specific_runs=[]):
    print 'Appending fixed probe runs.'
    print 'Appending run ' + str(all_runs[0]) + '.',
    single_runs = {}
    fp_moment_df = pd.read_hdf(file_path, key='run_'+str(all_runs[0])+'_moment_df')
    if all_runs[0] in specific_runs: single_runs[all_runs[0]] = fp_moment_df.copy()
    for run in all_runs[1:]:
        print '\rAppending run ' + str(run) + '.',
        try:
            temp_df = pd.read_hdf(file_path, key='run_'+str(run)+'_moment_df')
            fp_moment_df = fp_moment_df.append(temp_df)
            if run in specific_runs: single_runs[run] = temp_df.copy()
        except KeyError:
            print '\rNo run ' + str(run) + ' in hdf5 file.',
    print '\nDone appending fixed probe runs.'
    if specific_runs: return fp_moment_df, single_runs
    else: return fp_moment_df

def import_tr_moment_runs(file_path, all_runs, specific_runs=[], corrected=False):
    if corrected: print 'Appending corrected trolley runs.'
    else: print 'Appending un-corrected trolley runs.'
    print 'Appending run ' + str(all_runs[0]) + '.'
    single_runs = {}
    temp_df = pd.read_hdf(file_path, key='run_' + str(all_runs[0]) + '_moment_df')
    if corrected:
        tr_output_df = trfp.remove_trolley_effect(temp_df)
        if all_runs[0] in specific_runs: single_runs[all_runs[0]] = tr_output_df.copy()
    else:
        tr_output_df = temp_df.copy()
        if all_runs[0] in specific_runs: single_runs[all_runs[0]] = tr_output_df.copy()
    for run in all_runs[1:]:
        print 'Appending run ' + str(run) + '.'
        temp_df = pd.read_hdf(file_path, key='run_' + str(run) + '_moment_df')
        if corrected:
            temp_df = trfp.remove_trolley_effect(temp_df)
            if run in specific_runs: single_runs[run] = temp_df.copy()
            tr_output_df = tr_output_df.append(temp_df)
        else:
            if run in specific_runs: single_runs[run] = temp_df.copy()
            tr_output_df = tr_output_df.append(temp_df)
    if corrected: print '\nDone appending corrected trolley runs.'
    else: print '\nDone appending un-corrected trolley runs.'
    if specific_runs: return tr_output_df, single_runs
    else: return tr_output_df

# the NEW trolley station averaging routine

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
    tr_baseline = np.full((72,6), np.nan)
    fp_baseline = np.full((72,6), np.nan)
    summed_azimuth = np.full((72,6), np.nan)
    summed_pts = np.full((72,6), np.nan)
    baseline_time = np.full((72,6), np.nan)

    st6 = 0
    for st in range(72):
        
        num_probes = len(trfp.STATION_PROBE_ID[st])
        # first do m1-4 for all stations
        
        if station_edges[st+1] > station_edges[st]:
            mask = (corrected_df['tr_phi'] >= station_edges[st]) & (corrected_df['tr_phi'] < station_edges[st+1])
        else:  # case where we go over the 360 deg line
            mask = (corrected_df['tr_phi'] >= station_edges[st]) | (corrected_df['tr_phi'] < station_edges[st+1])

        out_df = corrected_df[mask].copy()
        summed_pts[st, :] = out_df.shape[0]
        summed_azimuth[st, :] = sum(out_df['tr_extent'].values)
        baseline_time[st, :] = sum(out_df.index.values)/summed_pts[st, 0]

        for m in range(6):
            st_id = 'tr,m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                tr_baseline[st, m] = sum(out_df['tr_extent'].values*out_df[st_id].values)/sum(out_df['tr_extent'].values)

            st_id = 'st'+str(st)+',m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                fp_baseline[st, m] = np.mean(out_df[st_id])
        
#         if num_probes == 4: continue  # moves to next iteration for 4 probe stations
        
#         # next do m5+ for all 6-probe stations

#         if station_edges_6_probe[st6+1] > station_edges_6_probe[st6]:
#             mask = (corrected_df['tr_phi'] >= station_edges_6_probe[st6]) & (corrected_df['tr_phi'] < station_edges_6_probe[st6+1])
#         else:  # case where we go over the 360 deg line
#             mask = (corrected_df['tr_phi'] >= station_edges_6_probe[st6]) | (corrected_df['tr_phi'] < station_edges_6_probe[st6+1])
        
#         out_df = corrected_df[mask].copy()
#         summed_pts[st, 4:6] = out_df.shape[0]
#         summed_azimuth[st, 4:6] = sum(out_df['tr_extent'].values)
#         baseline_time[st, 4:6] = sum(out_df.index.values)/summed_pts[st, 4]

#         for m in range(4,6):
#             st_id = 'tr,m'+str(m+1)
#             if sum(out_df['tr_extent'].values) != 0:
#                 tr_baseline[st, m] = sum(out_df['tr_extent'].values*out_df[st_id].values)/sum(out_df['tr_extent'].values)

#             st_id = 'st'+str(st)+',m'+str(m+1)
#             if sum(out_df['tr_extent'].values) != 0:
#                 fp_baseline[st, m] = np.mean(out_df[st_id])
        
#         st6 +=1
    
    return tr_baseline, fp_baseline, baseline_time, summed_azimuth, summed_pts

# generate virutal trolley measurements

def vtm_calc(fp_moment_df,
             baseline_time_1, baseline_time_2,
             tr_baseline_1, tr_baseline_2,
             fp_baseline_1, fp_baseline_2):
    
    vtm_df = fp_moment_df.copy()
#     vtm_baseline = fp_moment_df.copy()  # for old J method
    
    for st in range(72):
        num_probes = trfp.STATION_PROBE_NUMBER[st]
        
        if num_probes == 4:
            num_moments = 4
            if st == 41: J = trfp.J_4_PROBE_ST41
            elif st == 37 | st == 39: J = trfp.J_4_PROBE_ST37_ST39
            else: J = trfp.J_4_PROBE
            J = J[0:4, 0:4]
        else:
            num_moments = 5
            if st < 7:
                J = trfp.J_6_PROBE_OFFSET
            else: J = trfp.J_6_PROBE
        
        # first subtract fixed probe baseline from vtm_df
        for m in range(num_probes):
            stm = 'st'+str(st)+',m'+str(m+1)

            def __backwards_correction(time):
                c1 = fp_baseline_1[st, m]
                c2 = fp_baseline_2[st, m]
                t1 = baseline_time_1[st, m]
                t2 = baseline_time_2[st, m]
                return (c2-c1)/(t2-t1)*(time-t1) + c1

            vtm_df[stm] = vtm_df[stm] - __backwards_correction(vtm_df.index.values)
#             vtm_baseline[stm] = vtm_baseline[stm] - __backwards_correction(vtm_df.index.values)  # for old J method

        # next apply the Jacobian to the station
        stms = ['st'+str(st)+',m'+str(m+1) for m in range(num_moments)]
        vtm_df[stms] = vtm_df[stms].dot(np.transpose(J))
        
### The following is the OLD (non-vectorized) method for applying Jacobian.
### It requires a "baseline" dataframe that doesn't change values as the jacobian is applied.
#         for m in range(num_moments):
#             vtr_stm = 'st'+str(st)+',m'+str(m+1)
#             fp_stm = ['st'+str(st)+',m'+str(fp_m+1) for fp_m in np.arange(num_moments)]
#             vtm_df[vtr_stm] = vtm_baseline[fp_stm].dot(J[m])

        # finally add trolley baseline to vtm_df
        for m in range(num_probes):
            stm = 'st'+str(st)+',m'+str(m+1)
            
            def __backwards_correction(time):
                c1 = tr_baseline_1[st, m]
                c2 = tr_baseline_2[st, m]
                t1 = baseline_time_1[st, m]
                t2 = baseline_time_2[st, m]
                return (c2-c1)/(t2-t1)*(time-t1) + c1

            vtm_df[stm] = vtm_df[stm] + __backwards_correction(vtm_df.index.values)

    return vtm_df

# trolley footprint removal

def __split_by_nan(input_array):
    return [input_array[clump] for clump in np.ma.clump_unmasked(np.ma.masked_invalid(input_array))]

# def trolley_footprint_replacement(moment_df, veto_low=False, veto_high=False):

#     nomask_df = moment_df.copy()
#     mask_df = nomask_df.copy()
#     temp_avg_df = pd.DataFrame(index=mask_df.index)

#     split_station = []
#     all_good_stations = np.arange(6,72)  # not using the inflector stations
#     no_ground_loop_stations = np.array(range(6,16)+range(64,72))  # valid for 25 deg veto
    
#     if not veto_low:
#         barcode = trfp.STATION_BARCODE_PHI
#         veto_extent = 25
#         veto_low = (np.array(barcode)-(veto_extent-3)/2)%360
#         veto_high = (np.array(barcode)+3+(veto_extent-3)/2)%360

#     # first need to mask when trolley is near each station
#     for st in range(72):
#         stms = ['st' + str(st) + ',m' + str(m+1) for m in range(6)]

#         if veto_low[st] < veto_high[st]:
#             mask = (nomask_df['tr_phi']>veto_low[st]) & (nomask_df['tr_phi']<veto_high[st])
#         else:  # this happens when wrapping around 360 deg
#             mask = (nomask_df['tr_phi']>veto_low[st]) | (nomask_df['tr_phi']<veto_high[st])

#         if mask.iloc[0] & mask.iloc[-1]: split_station += [True]
#         else: split_station += [False]

#         mask_df[stms] = nomask_df[stms].mask(mask)

#         # next need to average all good stations that are not within 3 of current station
#         if st not in range(16, 23):  # note that these ranged were chosen for 25 deg veto
#             averaging_stations = np.delete(all_good_stations,
#                                            np.argwhere((np.abs((all_good_stations - st)%72)<=3)
#                                                        | (np.abs((all_good_stations - st)%72)>=69))
#                                           )
#         else:
#             averaging_stations = np.delete(no_ground_loop_stations,
#                                            np.argwhere((np.abs((no_ground_loop_stations - st)%72)<=3)
#                                                        | (np.abs((no_ground_loop_stations - st)%72)>=69))
#                                           )
#         for m in range(6):  # this will need to go over all moments
#             stm = 'st' + str(st) + ',m' + str(m+1)
#             avg_stms = ['st'+str(avg_st)+',m'+str(m+1) for avg_st in averaging_stations]
#             temp_avg_df[stm] = nomask_df[avg_stms].mean(axis=1).mask(~mask)

#     replaced_df = mask_df.copy()


#     # next need to remove the ring wide drift and replace with the station drift
#     for st in range(72):
#         num_moments = len(trfp.STATION_PROBE_ID[st])
#         for m in range(num_moments):
#             stm = 'st' + str(st) + ',m' + str(m+1)
#             num_endpts = 5
#             if not split_station[st]:

#                 inner_splits = __split_by_nan(temp_avg_df[stm].values)
#                 outer_splits = __split_by_nan(mask_df[stm].values)
#                 first_inner_avg = np.mean(inner_splits[0][0:num_endpts])
#                 last_inner_avg = np.mean(inner_splits[0][-num_endpts:])
#                 first_outer_avg = np.mean(outer_splits[0][-num_endpts:])
#                 last_outer_avg = np.mean(outer_splits[1][0:num_endpts])  # these all use 5 values to make extrapolation easier

#                 inner_delta_y = 0.5 * (last_inner_avg-first_inner_avg)/inner_splits[0].size * (num_endpts-1)
#                 outer_delta_y = 0.5 * (last_outer_avg-first_outer_avg)/inner_splits[0].size * (num_endpts-1)
#                 inner_lin_fit = np.linspace(first_inner_avg-inner_delta_y, last_inner_avg+inner_delta_y, num=inner_splits[0].size)
#                 outer_lin_fit = np.linspace(first_outer_avg+outer_delta_y, last_outer_avg-outer_delta_y, num=inner_splits[0].size)

#                 replacement_values = inner_splits[0] - inner_lin_fit + outer_lin_fit
#                 replaced_df[stm][replaced_df[stm].isna()] = np.array(replacement_values)

#             else:

#                 inner_splits = __split_by_nan(temp_avg_df[stm].values)
#                 first_inner_avg = [np.mean(inner_splits[0][0:num_endpts]), np.mean(inner_splits[1][0:num_endpts])]
#                 last_inner_avg = [np.mean(inner_splits[0][-num_endpts:]), np.mean(inner_splits[1][-num_endpts:])]
#                 inner_delta_y = [0.5 * (last_inner_avg[0]-first_inner_avg[0])/inner_splits[0].size * (num_endpts-1),
#                                  0.5 * (last_inner_avg[1]-first_inner_avg[1])/inner_splits[1].size * (num_endpts-1)]

#                 first_inner_fit = np.linspace(first_inner_avg[0]-inner_delta_y[0],
#                                               last_inner_avg[0]+inner_delta_y[0],
#                                               num=inner_splits[0].size)
#                 second_inner_fit = np.linspace(first_inner_avg[1]-inner_delta_y[1],
#                                               last_inner_avg[1]+inner_delta_y[1],
#                                               num=inner_splits[1].size)

#                 # use 260 seconds of data after (or before) the vetoed window to make a linear fit to approximate "station drift"
#                 outer_splits = __split_by_nan(mask_df[stm].values)
#                 dt = 1  # the time step, usually 1 sec, but might as well make it a variable
#                 num_pts = 260//dt  # integer period of 130 sec signal
#                 first_outer_fit_coeffs = np.polyfit(np.arange(num_pts)*dt, outer_splits[0][0:num_pts], deg=1)
#                 second_outer_fit_coeffs = np.polyfit(np.arange(num_pts)*dt, outer_splits[0][-num_pts:], deg=1)
#                 first_outer_fit = np.polyval(first_outer_fit_coeffs,
#                                              np.linspace(inner_splits[0].size*-dt, -dt, num=inner_splits[0].size)
#                                             )
#                 second_outer_fit = np.polyval(second_outer_fit_coeffs,
#                                               np.linspace(num_pts*dt, (num_pts+inner_splits[0].size)*dt,
#                                                           num=inner_splits[1].size)
#                                              )

#                 first_replacement_values = inner_splits[0] - first_inner_fit + first_outer_fit
#                 second_replacement_values = inner_splits[1] - second_inner_fit + second_outer_fit
#                 replacement_values = np.append(first_replacement_values, second_replacement_values)
#                 replaced_df[stm][replaced_df[stm].isna()] = np.array(replacement_values)

#     return replaced_df

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

        for m in range(6):

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
