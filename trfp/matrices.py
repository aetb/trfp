#!/usr/bin/env python
# coding: utf-8

# # Matrices
# 
# Contains definitions of the matrices for the matrix method

# In[ ]:


"""A module that contains the various matrices used in the matrix method."""

import numpy as np
import scipy.optimize
import trfp


# ## Helper functions
# Dunder function, can't (easily) be used outside of this script.
# Maybe replace with a helper function file?

# In[ ]:


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


# ## Moment matrices
# Transform vectors of frequencies into vectors of moments.

# In[ ]:


THETA_FP_6 = np.array([np.array([1, 1, 1, 1, 1, 1])/6.,  # dipole
                       np.array([1, 0, -1, 1, 0, -1])/-12.*4.5,  # n quad
                       np.array([1, 1, 1, -1, -1, -1])/46.2*4.5,  # s quad
                       np.array([1, 0, -1, -1, 0, 1])/-92.4/2*4.5**2,  # s sext
                       np.array([1, -2, 1, 1, -2, 1])/18./2*4.5**2,  # n sext
                       np.array([1, -2, 1, -1, 2, -1])/-138.6*4.5**3]  # NOT oct, no idea what this is
                     )

THETA_FP_4 = np.array([np.array([1, 0, 1, 0])/2.,  # dipole
                       np.array([1, -1, 1, -1])/-6.*4.5,  # n quad
                       np.array([1, 1, -1, -1])/30.8*4.5,  # s quad
                       np.array([1, -1, -1, 1])/-46.2/2*4.5**2]  # sext
                     )

# four probe stations in the garage (stations 37, 39, 41) have a different arrangement of probes
# STATION 37, 39
# TM, TO, BI, BM
THETA_FP_4_ST37_ST39 = np.array([np.array([1, 0, 0, 1])/2.,  #dipole
                                 np.array([-1, 1, -1, 1])/6.*4.5,  # n quad
                                 np.array([1, 1, -1, -1])/30.8*4.5,  # s quad
                                 np.array([-1, 1, 1, -1])/46.2/2*4.5**2]  # sext?
                               )
# STATION 41
# TO, TM, BI, BM
THETA_FP_4_ST41 = np.array([np.array([0, 1, 0, 1])/2.,  #dipole
                            np.array([1, -1, -1, 1])/6.*4.5,  # n quad
                            np.array([1, 1, -1, -1])/30.8*4.5,  # s quad
                            np.array([1, -1, 1, -1])/46.2/2*4.5**2]  # sext?
                          )

# No longer attempt to calculate either 18-pole.
# Jonathan may have changed this by accident (before 7/21/2020)

_MULTIPOLE_ORDER = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
_MULTIPOLE_SKEW = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1]
_MULTS = np.array([__multipole(_MULTIPOLE_ORDER[i], _MULTIPOLE_SKEW[i], 1, trfp.TR_X, trfp.TR_Y) for i in range(14)])
_MULTS[np.abs(_MULTS) < 1.0e-9] = 0 
THETA_TR = np.linalg.pinv(np.transpose(_MULTS))
THETA_TR = np.insert(THETA_TR, 12, np.zeros(17), 0)
THETA_TR = np.append(THETA_TR, np.zeros([2,17]), 0)
THETA_TR[np.abs(THETA_TR) < 1.0e-9] = 0


#_MULTIPOLE_ORDER = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
#_MULTIPOLE_SKEW = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1]
#_MULTS = np.array([__multipole(_MULTIPOLE_ORDER[i], _MULTIPOLE_SKEW[i], 1, trfp.TR_X, trfp.TR_Y) for i in range(14)])
#_MULTS[np.abs(_MULTS) < 1.0e-9] = 0 
#THETA_TR = np.linalg.pinv(np.transpose(_MULTS))
#THETA_TR = np.insert(THETA_TR, 12, np.zeros(16), 0)
#THETA_TR = np.append(THETA_TR, np.zeros([2,16]), 0)
#THETA_TR[np.abs(THETA_TR) < 1.0e-9] = 0
#theta = np.array([np.insert(THETA_TR[0], 8, 0)])

#for i in range(len(THETA_TR)-1):
#    row = np.array([np.insert(THETA_TR[i+1], 8, 0)]) #0 is added to the element indexed at8, 9th probe
#    theta = np.concatenate((theta, row), 0)
#THETA_TR = theta

# ## Jacobian matrices
# New Jacobians calculated analytically (../jacobian_analytic_calc.ipynb).
# 
# Note that stations 37, 39, and 41 have the same Jacobian because their geometry is the same. They are both recorded here just to emphasize that.

# In[ ]:


J_6_PROBE = np.array([[1.0, 0.0, 0.0, 0.0, 2.631605],
                      [0.0, 1.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 1.0]])

J_6_PROBE_OFFSET = np.array([[1.0, 0.222222, 0.0, 0.0, 2.680987],
                             [0.0, 1.0, 0.0, 0.0, 0.444444],
                             [0.0, 0.0, 1.0, 0.444444, 0.0],
                             [0.0, 0.0, 0.0, 1.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 1.0]])
#Uncomment if working with Station 5
# J_6_PROBE = np.array([[1.0, 0.0, 0.0, 0.0, 0.0],
#                       [0.0, 1.0, 0.0, 0.0, 0.0],
#                       [0.0, 0.0, 1.0, 0.0, 0.0],
#                       [0.0, 0.0, 0.0, 1.0, 0.0],
#                       [0.0, 0.0, 0.0, 0.0, 1.0]])

# J_6_PROBE_OFFSET = np.array([[1.0, 0.0, 0.0, 0.0, 0.0],
#                              [0.0, 1.0, 0.0, 0.0, 0.0],
#                              [0.0, 0.0, 1.0, 0.0, 0.0],
#                              [0.0, 0.0, 0.0, 1.0, 0.0],
#                              [0.0, 0.0, 0.0, 0.0, 1.0]])

J_4_PROBE = np.array([[1.0, 0.0, 0.0, 0.0, 2.927901],
                      [0.0, 1.0, 0.0, 0.0, -0.666667],
                      [0.0, 0.0, 1.0, -0.666667, 0.0],
                      [0.0, 0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 1.0]])

J_4_PROBE_ST41 = np.array([[1.0, 0.0, 0.0, 0.0, 2.927901],
                           [0.0, 1.0, 0.0, 0.0, 0.0],
                           [0.0, -0.194805, 1.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 1.0, -0.194805],
                           [0.0, 0.0, 0.0, 0.0, 1.0]])

J_4_PROBE_ST37_ST39 = np.array([[1.0, 0.0, 0.0, 0.0, 2.927901],
                                [0.0, 1.0, 0.0, 0.0, 0.0], 
                                [0.0, -0.194805, 1.0, 0.0, 0.0], 
                                [0.0, 0.0, 0.0, 1.0, -0.194805], 
                                [0.0, 0.0, 0.0, 0.0, 1.0]])

