#This does only trolley probe dropping, Analysis Helper2 does fixed probe dropping
import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz


import gm2
import trfp
import helper_function_candidates as helper


def __multipole(order, skew, strength, x_pos, y_pos):
    """Returns the magnitude of a B-field given multipole parameters
    and position.
    """
    # Takes a multipole strength A normalized to 4.5 cm
    r_pos = np.sqrt(x_pos**2 + y_pos**2)
    theta = np.arctan2(y_pos, x_pos)
    if skew == 0:
        b_magnitude = strength * (r_pos/4.5)**order * np.cos(order*theta)
    if skew == 1:
        b_magnitude = strength * (r_pos/4.5)**order * np.sin(order*theta)
    return b_magnitude
# Probe Geometry.
class Geometry:
    def __init__(self, x, y, p):
        self.xpos = x
        self.ypos = y
        self.probes = p.astype('int')
        if len(self.xpos) != len(self.ypos):
            print("In Geometry length of xpos and ypos are not equal")
        elif len(self.ypos) != len(self.probes):
            print("In Geometry length of probes and ypos are not equal")
    def IsReduced(self):
        isReduced = True
        for i in self.probes:
            if i == 0:
                isReduced = i
        return isReduced
    
def DropPos(geom):
    badProbes = np.array([])
    for i in range(len(geom.probes)):
        if not geom.probes[i]:
            badProbes = np.append(badProbes, i)
    badProbes = badProbes.astype(int)
    X = geom.xpos.copy()
    Y = geom.ypos.copy()
    px_new = np.delete(X, badProbes) # Only the x pos of working probes stay
    py_new = np.delete(Y, badProbes) # Same for y
    PROBES_NEW = np.delete(geom.probes, badProbes) #Now this is an array of 1s.
    reducedGeometry = Geometry(px_new, py_new, PROBES_NEW)
    return reducedGeometry, badProbes

def Theta(geom, badProbes, mult_order, mult_skew):    
    #"If then" for truncating moments, can recover all moments from dropping 1 trfp 
    if not geom.IsReduced():
        print("The geometry given has non-working probes")
        return
    lenOrigProbes = len(geom.probes) + len(badProbes)
    if len(badProbes) + len(mult_order) > lenOrigProbes:
        _MULTS_FP_N = np.array([__multipole(mult_order[i], mult_skew[i], 1, geom.xpos, geom.ypos)\
                                for i in range(len(mult_order) - len(badProbes))])
    else:
        _MULTS_FP_N = np.array([__multipole(mult_order[i], mult_skew[i], 1, geom.xpos, geom.ypos)\
                                for i in range(len(mult_order))])
        
    _MULTS_FP_N[np.abs(_MULTS_FP_N) < 1.0e-9] = 0
    
    INV_MULT = np.linalg.pinv(np.transpose(_MULTS_FP_N))
    INV_MULT[np.abs(INV_MULT) < 1.0e-9] = 0
    
    I = np.dot(INV_MULT,np.transpose(_MULTS_FP_N))
    I[np.abs(I) < 1.0e-9] =0
    #print(I)
    #Add the columns of 0
    inv_mult = np.zeros((len(INV_MULT), lenOrigProbes))
    j = 0
    if len(badProbes) > 0:
        badProbesCopy = badProbes.copy()
        for i in range(lenOrigProbes):
            if i != badProbesCopy[0]:
                inv_mult[:,i] = INV_MULT[:,j]
                j+=1
            elif (i == badProbesCopy[0]) and (len(badProbesCopy) > 1):
                badProbesCopy = np.delete(badProbesCopy, 0)
                  
        #Add rows at the end, only the first moments can be calculated
        for i in range(lenOrigProbes - len(geom.probes)):
            inv_mult = np.concatenate((inv_mult, np.zeros((1,lenOrigProbes))))
        return(inv_mult)
    else:
        return(INV_MULT)

    
def ProbeDrops(geom):
    reducedGeometry, badProbes = DropPos(geom)# Return another geom
    
    _MULTIPOLE_ORDER = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
    _MULTIPOLE_SKEW = [0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1]
    #Check if station or Trolley
    Trolley = False
    if (len(geom.probes)) > 6:
        Trolley = True
    
    if Trolley:
        _MULTIPOLE_SKEW = _MULTIPOLE_SKEW[:len(geom.probes)]
        _MULTIPOLE_ORDER = _MULTIPOLE_ORDER[:len(geom.probes)]
    elif len(badProbes) > 1:
        _MULTIPOLE_SKEW = _MULTIPOLE_SKEW[:len(geom.probes)]
        _MULTIPOLE_ORDER = _MULTIPOLE_ORDER[:len(geom.probes)]
    else:
        _MULTIPOLE_SKEW = _MULTIPOLE_SKEW[:5] #Up to m5
        _MULTIPOLE_ORDER = _MULTIPOLE_ORDER[:5] #Up to m5

    return Theta(reducedGeometry, badProbes, _MULTIPOLE_ORDER, _MULTIPOLE_SKEW)
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
    

# ## Root to interpolated data frame
# `root_to_pandas`

# In[ ]:


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


def calc_moment_df(interp_df, trolleyProbes = np.zeros(17)+1,\
                   fixedProbes = np.zeros((76,6))+1):
    #Make regular regular 4 and 6 FP thetas -JS
    FP6Geometry = Geometry(trfp.FP6_X, trfp.FP6_Y, np.zeros(6)+1)
    FP6Theta = ProbeDrops(FP6Geometry)
        
    FP4Geometry = Geometry(trfp.FP4_X, trfp.FP4_Y, np.zeros(4)+1)
    FP4Theta = ProbeDrops(FP4Geometry)
    
    tr_run = 'tr_phi' in interp_df.columns
    
    moment_df = pd.DataFrame(index=interp_df.index)
    
    # calculate trolley moments if needed
    if tr_run:
        #Trolley Probe Dropping -JS
        TrolleyGeometry = Geometry(trfp.TR_X, trfp.TR_Y, trolleyProbes)
        TrolleyTheta = ProbeDrops(TrolleyGeometry)
        TrolleyTheta = np.insert(TrolleyTheta, 12, np.zeros(17), axis = 0)
        TrolleyTheta = np.concatenate((TrolleyTheta, np.zeros((2,17))))
        
        moment_df['tr_phi'] = interp_df['tr_phi'].copy()
        print 'Calculating trolley moments.'
        theta_tr = TrolleyTheta
        for m in np.arange(17):
            tr_probes = ['tr'+str(probe) for probe in np.arange(17)]
            moment_df['tr,m'+str(m+1)] = interp_df[tr_probes].dot(theta_tr[m])

    # create the 72*6 fixed probe moments
    for station in np.arange(72):
        print '\rCalculating station ' + str(station) + ' moments.',
        fp_st = ['fp'+str(fp) for fp in trfp.STATION_PROBE_ID[station]]
        
        # Check number of probes in station and rearrange geometries -JS
        if len(trfp.STATION_PROBE_ID[station]) == 6:
            if not( 0 in fixedProbes[station]):
                theta_fp = FP6Theta
                #theta_fp = trfp.THETA_FP_6
            else: 
                StGeom = Geometry(trfp.FP6_X, trfp.FP6_Y, fixedProbes[station])
                StTheta = ProbeDrops(StGeom)
                theta_fp = StTheta
        
        elif(station == 41):
            St41Geom = Geometry(trfp.FP4_X_ST41, trfp.FP4_Y, fixedProbes[station][:4])
            St41Theta = ProbeDrops(St41Geom)
            theta_fp = St41Theta
        
        elif (station == 37) or (station == 39):
            St37_39Geom = Geometry(trfp.FP4_X_ST37_ST39, trfp.FP4_Y,\
                                   fixedProbes[station][:4])
            St37_39Theta = ProbeDrops(St37_39Geom)
            theta_fp = St37_39Theta
            
        else:
            if not (0 in fixedProbes[station][:4]):
                theta_fp = FP4Theta
            else:
                StGeom = Geometry(trfp.FP4_X, trfp.FP4_Y, fixedProbes[station][:4])
                StTheta = ProbeDrops(StGeom)
                theta_fp = StTheta

        # # choose proper theta matrix     #Not needed when probe dropping -JS
        # theta_fp = _choose_theta(station)#

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

# generate virutal trolley measurements

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