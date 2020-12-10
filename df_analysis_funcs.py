import numpy as np
import pandas as pd

import analysis_helper as helper
import helper_function_candidates as helper_old

import gm2
import trfp


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

def _get_keys(filename):
    hdf = pd.HDFStore(filename, mode='r')
    keys = hdf.keys()
    hdf.close()
    output_keys = []
    for key in keys:
        if key != '/subrun_df': output_keys.append(key[1:])
    return output_keys

def read_dfs(filename):
    keys = _get_keys(filename)
    interp_dfs = {}
    for key in keys:
        print '\nReading ' + key
        interp_dfs[key] = pd.read_hdf(filename, key=key)
    subrun_df = pd.read_hdf(filename, key='subrun_df')
    return interp_dfs, keys, subrun_df

def interp_to_moment(interp_dfs, keys):
    moment_dfs = {}
    for key in keys:
        print '\nCalculating moments for ' + key
        moment_dfs[key] = helper.calc_moment_df(interp_dfs[key])
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
            corrected_dfs[key] = helper_old.trolley_footprint_replacement(moment_dfs[key])
        elif key[:2] == 'fp': corrected_dfs[key] = moment_dfs[key].copy()
    print '\n'
    return corrected_dfs

def bloch_style_moments(corrected_dfs, keys):
    print '\nImplementing Bloch-style treatment of stations 1, 3, 5, and 54'
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
            tr_baselines[key], fp_baselines[key], baseline_times[key], summed_azimuths[key], summed_pts[key] = helper_old.trolley_run_station_average(corrected_dfs[key])

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
            vtm_dfs[key] = helper.vtm_calc(corrected_dfs[key],
                                           baselines['time'][pair_dict[key][0]], baselines['time'][pair_dict[key][1]],
                                           baselines['tr'][pair_dict[key][0]], baselines['tr'][pair_dict[key][1]],
                                           baselines['fp'][pair_dict[key][0]], baselines['fp'][pair_dict[key][1]]
                                           )
    
    return vtm_dfs
