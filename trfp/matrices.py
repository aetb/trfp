"""A module that contains classes which generate
the various matrices used in the matrix method.
"""

import numpy as np
import trfp


class FieldToMomentMatrix(object):

    """A class that generate matrices taking raw probe
    measurements to field moments.
    """

    def __init__(self):
        self.theta_fp_6 = np.array(
            [np.array([1, 1, 1, 1, 1, 1])/6.,  # dipole
             np.array([1, 0, -1, 1, 0, -1])/-12.*4.5,  # n quad
             np.array([1, 1, 1, -1, -1, -1])/-46.2*4.5,  # s quad
             np.array([1, 0, -1, -1, 0, 1])/92.4*4.5**2,  # s sext
             np.array([1, -2, 1, 1, -2, 1])/18.*4.5**2,  # n sext
             np.array([1, -2, 1, -1, 2, -1])/-138.6*4.5**3]  # oct
        )

        self.theta_fp_4 = np.array(
            [np.array([1, 0, 1, 0])/2.,  # dipole
             np.array([1, -1, 1, -1])/-6.*4.5,  # n quad
             np.array([1, 1, -1, -1])/-30.8*4.5,  # s quad
             np.array([1, -1, -1, 1])/46.2*4.5**2]  # sext
        )

        self.order = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]
        self.skew = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]

        self.theta_tr = self.__theta_tr()

    def __multipole(self, order, skew, strength, x_pos, y_pos):

        # Takes a multipole strength A normalized to 4.5 cm
        r_pos = np.sqrt(x_pos**2 + y_pos**2)
        theta = np.arctan2(y_pos, x_pos)
        if skew == 0:
            b_magnitude = strength * (r_pos/4.5)**order * np.cos(order*theta)
        if skew == 1:
            b_magnitude = strength * (r_pos/4.5)**order * np.sin(order*theta)
        return b_magnitude

    def __theta_tr(self):
        geometry = trfp.Geometry()

        out = np.transpose(np.array(
            [self.__multipole(self.order[i], self.skew[i], 1, geometry.tr_x, geometry.tr_y)
             for i in np.arange(17)]))
        theta_tr = np.linalg.pinv(out)
        return theta_tr
