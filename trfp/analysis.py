"""Contains functions for analyzing trolley and fixed probe runs."""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import seaborn as sns

import gm2
import trfp

from datetime import datetime
import pytz

def remove_trolley_effect(trolley_run):
    '''DOC STRING'''
    barcode = trfp.STATION_BARCODE_PHI
    
    trolley_effect_removed_df = trolley_run.moment_df.copy()
    
    for st in range(72):
        print '\rRemoving trolley image from station ' + str(st) + '.',
        for m in range(1,7):
            st_m = 'st' + str(st) + ",m" + str(m)
            
            # Unwrap the fixed probe data versus trolley position
            raw_data = trolley_run.moment_df[['tr_phi', st_m]].copy()
            raw_low = raw_data.copy()
            raw_high = raw_data.copy()
            raw_low['tr_phi'] = raw_low['tr_phi'] - 360
            raw_high['tr_phi'] = raw_high['tr_phi'] + 360
            unwrap_nomask_df = pd.concat([raw_low, raw_data, raw_high])
            
            unwrap_mask_df = unwrap_nomask_df.copy()
            mask = ((unwrap_nomask_df['tr_phi']>barcode[st]-2) & (unwrap_nomask_df['tr_phi']<barcode[st]+5) |
                    (unwrap_nomask_df['tr_phi']>barcode[st]-2-360) & (unwrap_nomask_df['tr_phi']<barcode[st]+5-360) |
                    (unwrap_nomask_df['tr_phi']>barcode[st]-2+360) & (unwrap_nomask_df['tr_phi']<barcode[st]+5+360))
            unwrap_mask_df[st_m] = unwrap_nomask_df[st_m].mask(mask)
            unwrap_mask_df['tr_phi'] = unwrap_nomask_df['tr_phi']
            
            unwrap_filled_df = unwrap_mask_df.copy()
            temp = unwrap_filled_df.rolling(int(100),win_type='triang',min_periods=1,center=True).mean()
            unwrap_filled_df[st_m] = unwrap_filled_df[st_m].mask(mask, temp[st_m])
            
            length = raw_data.shape[0]
            filled_df = unwrap_filled_df.iloc[length:2*length,:]
            
            trolley_effect_removed_df[st_m] = filled_df[st_m]
    
    print '\rFinished removing trolley images from ' + str(length) + ' events.'
    return trolley_effect_removed_df