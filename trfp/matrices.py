"""A module that contains the various matrices used in the matrix method."""

import numpy as np
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

THETA_FP_6 = np.array([np.array([1, 1, 1, 1, 1, 1])/6.,  # dipole
                       np.array([1, 0, -1, 1, 0, -1])/-12.*4.5,  # n quad
                       np.array([1, 1, 1, -1, -1, -1])/-46.2*4.5,  # s quad
                       np.array([1, 0, -1, -1, 0, 1])/92.4*4.5**2,  # s sext
                       np.array([1, -2, 1, 1, -2, 1])/18.*4.5**2,  # n sext
                       np.array([1, -2, 1, -1, 2, -1])/-138.6*4.5**3]  # oct
                     )

THETA_FP_4 = np.array([np.array([1, 0, 1, 0])/2.,  # dipole
                       np.array([1, -1, 1, -1])/-6.*4.5,  # n quad
                       np.array([1, 1, -1, -1])/-30.8*4.5,  # s quad
                       np.array([1, -1, -1, 1])/46.2*4.5**2]  # sext
                     )

__MULTIPOLE_ORDER = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]
__MULTIPOLE_SKEW = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]

THETA_TR = np.linalg.pinv(np.transpose(np.array([__multipole(__MULTIPOLE_ORDER[i],
                                                             __MULTIPOLE_SKEW[i],
                                                             1, trfp.TR_X,
                                                             trfp.TR_Y
                                                            )
                                                 for i in np.arange(17)])))
