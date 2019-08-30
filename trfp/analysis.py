#!/usr/bin/env python
# coding: utf-8

# # Analysis
# Contains helper functions and main analysis functions.

# In[ ]:


"""Contains functions for analyzing trolley and fixed probe runs."""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import seaborn as sns

import gm2
import trfp


# ## Remove trolley footprint
# Function that removes the trolley footprint from fixed probe measurements.
# Uses ring-wide drift as relacement.

# In[ ]:


def remove_trolley_effect(trolley_moment_df):
    '''DOC STRING'''
    veto_extent = 25
    barcode = trfp.STATION_BARCODE_PHI
    
    trolley_effect_removed_df = trolley_moment_df.copy()
    
    for st in range(72):
        print '\rRemoving trolley image from station ' + str(st) + '.',
        for m in range(1,7):
            st_m = 'st' + str(st) + ",m" + str(m)
            
            # Unwrap the fixed probe data versus trolley position
            raw_data = trolley_moment_df[['tr_phi', st_m]].copy()
            raw_low = raw_data.copy()
            raw_high = raw_data.copy()
            raw_low['tr_phi'] = raw_low['tr_phi'] - 360
            raw_high['tr_phi'] = raw_high['tr_phi'] + 360
            unwrap_nomask_df = pd.concat([raw_low, raw_data, raw_high])
            
            unwrap_mask_df = unwrap_nomask_df.copy()
#             mask = ((unwrap_nomask_df['tr_phi']>barcode[st]-2) & (unwrap_nomask_df['tr_phi']<barcode[st]+5) |
#                     (unwrap_nomask_df['tr_phi']>barcode[st]-2-360) & (unwrap_nomask_df['tr_phi']<barcode[st]+5-360) |
#                     (unwrap_nomask_df['tr_phi']>barcode[st]-2+360) & (unwrap_nomask_df['tr_phi']<barcode[st]+5+360))
            veto_adjust = (veto_extent-7)/2
            mask = (  (unwrap_nomask_df['tr_phi']>barcode[st]-2-veto_adjust)
                       & (unwrap_nomask_df['tr_phi']<barcode[st]+5+veto_adjust)
                    | (unwrap_nomask_df['tr_phi']>barcode[st]-2-veto_adjust-360)
                       & (unwrap_nomask_df['tr_phi']<barcode[st]+5+veto_adjust-360)
                    | (unwrap_nomask_df['tr_phi']>barcode[st]-2-veto_adjust+360)
                       & (unwrap_nomask_df['tr_phi']<barcode[st]+5+veto_adjust+360))
            
            unwrap_mask_df[st_m] = unwrap_nomask_df[st_m].mask(mask)
            unwrap_mask_df['tr_phi'] = unwrap_nomask_df['tr_phi']
            
            unwrap_filled_df = unwrap_mask_df.copy()
            temp = unwrap_filled_df.rolling(int(500),win_type='triang',min_periods=1,center=True).mean()
            temp = temp.rolling(int(500),win_type='triang',min_periods=1,center=True).mean()
            unwrap_filled_df[st_m] = unwrap_filled_df[st_m].mask(mask, temp[st_m])
            
            length = raw_data.shape[0]
            filled_df = unwrap_filled_df.iloc[length:2*length,:]
            
            trolley_effect_removed_df[st_m] = filled_df[st_m]
    
    print '\rFinished removing trolley images from ' + str(length) + ' events.'
    return trolley_effect_removed_df


# ## Average fields station-wise during trolley runs
# Calculates average moments during trolley from for synchronization.

# In[ ]:


def trolley_run_station_average(corrected_df):
    station_phi = trfp.STATION_BARCODE_PHI
    station_edges = trfp.STATION_BARCODE_EDGES

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
    tr_baseline = np.empty([72,17])
    fp_baseline = np.empty([72,6])
    summed_azimuth = np.empty(72)
    summed_pts = np.empty(72)
    baseline_time = np.empty(72)
    tr_baseline[:] = np.nan
    fp_baseline[:] = np.nan
    summed_azimuth[:] = np.nan
    summed_pts[:] = np.nan
    baseline_time[:] = np.nan

    for st in range(72): 
        if station_edges[st+1] > station_edges[st]:
            mask = (corrected_df['tr_phi'] >= station_edges[st]) & (corrected_df['tr_phi'] < station_edges[st+1])
        else:  # case where we go over the 360 deg line
            mask = (corrected_df['tr_phi'] >= station_edges[st]) | (corrected_df['tr_phi'] < station_edges[st+1])

        out_df = corrected_df[mask]
        summed_pts[st] = out_df.shape[0]
        summed_azimuth[st] = sum(out_df['tr_extent'].values)        
        baseline_time[st] = sum(out_df.index.values)/summed_pts[st]

        for m in range(17):
            st_id = 'tr,m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                tr_baseline[st, m] = sum(out_df['tr_extent'].values*out_df[st_id].values)/sum(out_df['tr_extent'].values)
            else:
                tr_baseline[st, m] = np.nan
        for m in range(6):
            st_id = 'st'+str(st)+',m'+str(m+1)
            if sum(out_df['tr_extent'].values) != 0:
                fp_baseline[st, m] = np.mean(out_df[st_id])
            else:
                fp_baseline[st, m] = np.nan
    
    return tr_baseline, fp_baseline, baseline_time, summed_azimuth, summed_pts

