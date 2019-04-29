"""A module that contains the various matrices used in the matrix method."""

import numpy as np
import scipy.optimize
import trfp

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

def __lin_fit(x, a, b):
    return a + b*x

THETA_FP_6 = np.array([np.array([1, 1, 1, 1, 1, 1])/6.,  # dipole
                       np.array([1, 0, -1, 1, 0, -1])/-12.*4.5,  # n quad
                       np.array([1, 1, 1, -1, -1, -1])/46.2*4.5,  # s quad
                       np.array([1, 0, -1, -1, 0, 1])/-92.4*4.5**2,  # s sext
                       np.array([1, -2, 1, 1, -2, 1])/18.*4.5**2,  # n sext
                       np.array([1, -2, 1, -1, 2, -1])/-138.6*4.5**3]  # NOT oct, no idea what this is
                     )

THETA_FP_4 = np.array([np.array([1, 0, 1, 0])/2.,  # dipole
                       np.array([1, -1, -1, 1])/-6.*4.5,  # n quad
                       np.array([1, 1, -1, -1])/30.8*4.5,  # s quad
                       np.array([1, -1, 1, -1])/46.2*4.5**2]  # sext
                     )

# four probe stations in the garage (Yoke , stations 37, 39, 41) have a different arrangement of probes
# STATION 37, 39
# TM, TO, BI, BM
THETA_FP_4_ST37_ST39 = np.array([np.array([1, 0, 0, 1])/2.,  #dipole
                              np.array([-1, 1, -1, 1])/6.*4.5,  # n quad
                              np.array([1, 1, -1, -1])/30.8*4.5,  # s quad
                              np.array([-1, 1, 1, -1])/46.2*4.5**2]  # sext?
                            )
# STATION 41
# TO, TM, BI, BM
THETA_FP_4_ST41 = np.array([np.array([0, 1, 0, 1])/2.,  #dipole
                              np.array([1, -1, -1, 1])/6.*4.5,  # n quad
                              np.array([1, 1, -1, -1])/30.8*4.5,  # s quad
                              np.array([1, -1, 1, -1])/46.2*4.5**2]  # sext?
                            )

__MULTIPOLE_ORDER = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]
__MULTIPOLE_SKEW = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]

THETA_TR = np.linalg.pinv(np.transpose(np.array([__multipole(__MULTIPOLE_ORDER[i],
                                                             __MULTIPOLE_SKEW[i],
                                                             1, trfp.TR_X,
                                                             trfp.TR_Y
                                                            )
                                                 for i in np.arange(17)])))

# the following are the Jacobians that take fixed probe to trolley

def __jacobian_calc(probes, offset=False, weird_station=-1):
    tr_x = trfp.TR_X
    tr_y = trfp.TR_Y
    if probes == 6:
        THETA_FP = THETA_FP_6
        if offset:
            fp_x = trfp.FP6_X_OFFSET
        else:
            fp_x = trfp.FP6_X
        fp_y = trfp.FP6_Y
    else:
        probes = 4
        if weird_station == 41:
            THETA_FP = THETA_FP_4_ST41
            fp_x = trfp.FP4_X_ST41
        elif (weird_station == 37) | (weird_station == 39):
            THETA_FP = THETA_FP_4_ST37_ST39
            fp_x = trfp.FP4_X_ST37_ST39
        else:
            THETA_FP = THETA_FP_4
            fp_x = trfp.FP4_X
        fp_y = trfp.FP4_Y
        
    As = np.arange(-10,10)

    dfp_dtr = np.zeros((probes, probes))
    
    for ii in np.arange(probes):
        N = __MULTIPOLE_ORDER[ii]
        s = __MULTIPOLE_SKEW[ii]

        tr_out = np.empty((len(As),probes))
        tr_out[:] = np.nan
        fp_out = np.empty((len(As),probes))
        fp_out[:] = np.nan

        for jj in np.arange(len(As)):
            A = As[jj]
            B_tr = __multipole(N, s, A, tr_x, tr_y)
            B_fp = __multipole(N, s, A, fp_x, fp_y)

            tr_out[jj,:] = np.matmul(THETA_TR, B_tr)[0:probes]
            fp_out[jj,:] = np.matmul(THETA_FP, B_fp)

        for kk in np.arange(probes):
            coeffs, covar = scipy.optimize.curve_fit(__lin_fit, As, fp_out[:,kk])
            dfp_dtr[kk,ii] = coeffs[1]
    return dfp_dtr

# the following are dtr/dfp

J_6_PROBE = np.linalg.inv(__jacobian_calc(6, False)[0:5, 0:5])

J_6_PROBE_OFFSET = np.linalg.inv(__jacobian_calc(6, True)[0:5, 0:5])

J_4_PROBE = np.linalg.inv(__jacobian_calc(4, False))

J_4_PROBE_ST37_ST39 = np.linalg.inv(__jacobian_calc(4, False, weird_station=37))
# These last 2 should be equal, because the geometry is the same.
J_4_PROBE_ST41 = np.linalg.inv(__jacobian_calc(4, False, weird_station=41))


