#!/usr/bin/env python
# coding: utf-8

# # Geometry
# Contains constants encoding the geometry of the g-2 probes and magnet.

# In[ ]:


"""This module constains constants about the geometry of the experiment."""

import numpy as np


# ## Probe positions
# In-slice positions of the fixed probes and trolley probes.

# In[ ]:


FP4_X = np.array([0, 3, 0, 3])

# STATION 41
# TO, TM, BI, BM
FP4_X_ST41 = np.array([3, 0, -3, 0])

# STATION 37, 39
# TM, TO, BI, BM
FP4_X_ST37_ST39 = np.array([0, 3, -3, 0])

FP4_Y = np.array([7.7, 7.7, -7.7, -7.7])

FP6_X = np.array([-3, 0, 3, -3, 0, 3])
FP6_X_OFFSET = np.array([-4, -1, 2, -4, -1, 2])
FP6_Y = np.array([7.7, 7.7, 7.7, -7.7, -7.7, -7.7])

TR_X = np.array([0]
                + [1.75 * np.sin(2*np.pi/4*i)
                   for i in np.arange(4)]
                + [3.5 * np.sin(2*np.pi/12*i)
                   for i in np.arange(12)])
TR_Y = np.array([0]
                + [-1.75 * np.cos(2*np.pi/4*i)
                   for i in np.arange(4)]
                + [-3.5 * np.cos(2*np.pi/12*i)
                   for i in np.arange(12)])


# ## Station-probe information
# Which probes are in which station.
# 
# Also `STATION_RIGHT_PHI`, which is the (incorrect) nominal positions of the stations.
# Don't use `STATION_RIGHT_PHI`. It's just here for posterity.

# In[ ]:


STATION_RING_PHI = [348, 352, 358, 3, 8, 12, 18, 22, 28,
                    33, 38, 42, 48, 52, 58, 63, 68, 72,
                    78, 82, 88, 93, 98, 102, 108, 112,
                    118, 123, 128, 132, 138, 142, 148,
                    153, 158, 162, 168, 172, 178, 183,
                    188, 192, 198, 202, 208, 213, 218,
                    222, 228, 232, 238, 243, 248, 252,
                    258, 262, 268, 273, 278, 282, 288,
                    292, 298, 303, 308, 312, 318, 322,
                    328, 333, 338, 342]

STATION_PROBE_ID = [[0, 1, 2, 18, 19, 20],  # 0
                    [3, 4, 5, 21, 22, 23],  # 1
                    [6, 7, 8, 24, 25, 26],  # 2
                    [9, 10, 11, 27, 28, 29],  # 3
                    [12, 13, 14, 30, 31, 32],  # 4
                    [15, 16, 17, 33, 34, 35],  # 5
                    [36, 37, 38, 51, 52, 53],  # 6
                    [39, 40, 54, 55],  # 7
                    [41, 42, 43, 56, 57, 58],  # 8
                    [44, 45, 59, 60],  # 9
                    [46, 47, 48, 61, 62, 63],  # 10
                    [49, 50, 64, 65],  # 11
                    [66, 67, 68, 81, 82, 83],  # 12
                    [69, 70, 84, 85],  # 13
                    [71, 72, 73, 86, 87, 88],  # 14
                    [74, 75, 89, 90],  # 15
                    [76, 77, 78, 91, 92, 93],  # 16
                    [79, 80, 94, 95],  # 17
                    [96, 97, 98, 114, 115, 116],  # 18
                    [99, 100, 101, 117, 118, 119],  # 19
                    [102, 103, 104, 120, 121, 122],  # 20
                    [105, 106, 107, 123, 124, 125],  # 21
                    [108, 109, 110, 126, 127, 128],  # 22
                    [111, 112, 113, 129, 130, 131],  # 23
                    [132, 133, 134, 150, 151, 152],  # 24
                    [135, 136, 137, 153, 154, 155],  # 25
                    [138, 139, 140, 156, 157, 158],  # 26
                    [141, 142, 143, 159, 160, 161],  # 27
                    [144, 145, 146, 162, 163, 164],  # 28
                    [147, 148, 149, 165, 166, 167],  # 29
                    [168, 169, 170, 183, 184, 185],  # 30
                    [171, 172, 186, 187],  # 31
                    [173, 174, 175, 188, 189, 190],  # 32
                    [176, 177, 191, 192],  # 33
                    [178, 179, 180, 193, 194, 195],  # 34
                    [181, 182, 196, 197],  # 35
                    [198, 199, 200, 213, 214, 215],  # 36
                    [201, 202, 216, 217],  # 37
                    [203, 204, 205, 218, 219, 220],  # 38
                    [206, 207, 221, 222],  # 39
                    [208, 209, 210, 223, 224, 225],  # 40
                    [211, 212, 226, 227],  # 41
                    [228, 229, 230, 243, 244, 245],  # 42
                    [231, 232, 246, 247],  # 43
                    [233, 234, 235, 248, 249, 250],  # 44
                    [236, 237, 251, 252],  # 45
                    [238, 239, 240, 253, 254, 255],  # 46
                    [241, 242, 256, 257],  # 47
                    [258, 259, 260, 273, 274, 275],  # 48
                    [261, 262, 276, 277],  # 49
                    [263, 264, 265, 278, 279, 280],  # 50
                    [266, 267, 281, 282],  # 51
                    [268, 269, 270, 283, 284, 285],  # 52
                    [271, 272, 286, 287],  # 53
                    [288, 289, 290, 303, 304, 305],  # 54
                    [291, 292, 306, 307],  # 55
                    [293, 294, 295, 308, 309, 310],  # 56
                    [296, 297, 311, 312],  # 57
                    [298, 299, 300, 313, 314, 315],  # 58
                    [301, 302, 316, 317],  # 59
                    [318, 319, 320, 333, 334, 335],  # 60
                    [321, 322, 336, 337],  # 61
                    [323, 324, 325, 338, 339, 340],  # 62
                    [326, 327, 341, 342],  # 63
                    [328, 329, 330, 343, 344, 345],  # 64
                    [331, 332, 346, 347],  # 65
                    [348, 349, 350, 363, 364, 365],  # 66
                    [351, 352, 366, 367],  # 67
                    [353, 354, 355, 368, 369, 370],  # 68
                    [356, 357, 371, 372],  # 69
                    [358, 359, 360, 373, 374, 375],  # 70
                    [361, 362, 376, 377]]  # 71

STATION_PROBE_NUMBER = [len(probes) for probes in STATION_PROBE_ID]


# ## Station barcode positions
# Position of the fixed probe stations.
# Calculated by using the trolley footprint.
# Lined up with the trolley probe active regions in the barcode basis.

# In[ ]:


STATION_BARCODE_PHI = [350.17, 354.33, 358.84, 4.34, 9.33,
                       14.23, 20.2, 23.32, 29.33, 34.33, 39.33,
                       43.31, 50.2, 53.31, 59.34, 64.33, 69.34,
                       73.3, 80.17, 83.32, 89.31, 94.33, 99.31,
                       103.3, 110.19, 113.31, 119.34, 124.32,
                       129.3, 133.27, 140.21, 143.38, 149.34,
                       154.33, 159.36, 163.31, 170.19, 173.39,
                       179.36, 184.35, 189.36, 193.31, 200.17,
                       203.38, 209.34, 214.31, 219.34, 223.32,
                       230.22, 233.42, 239.36, 244.34, 249.38,
                       253.34, 260.2, 263.47, 269.35, 274.34,
                       279.38, 283.41, 290.19, 293.37, 299.34,
                       304.31, 309.34, 313.29, 320.19, 323.41,
                       329.33, 334.35, 339.35, 343.41]

STATION_BARCODE_EDGES = (STATION_BARCODE_PHI+np.roll(STATION_BARCODE_PHI,1))/2
if STATION_BARCODE_EDGES[3] >= 180.:  # accounts for wrap around at station 3
    STATION_BARCODE_EDGES[3] = STATION_BARCODE_EDGES[3]-180.
else:
    STATION_BARCODE_EDGES[3] = STATION_BARCODE_EDGES[3]+180.
STATION_BARCODE_EDGES = np.append(STATION_BARCODE_EDGES, STATION_BARCODE_EDGES[0])

STATION_LIST_6_PROBES = []
for ii in range(72):
    if STATION_PROBE_NUMBER[ii] == 6:
        STATION_LIST_6_PROBES.append(ii)

STATION_BARCODE_PHI_6 = [STATION_BARCODE_PHI[st] for st in STATION_LIST_6_PROBES]

STATION_BARCODE_EDGES_6 = (STATION_BARCODE_PHI_6+np.roll(STATION_BARCODE_PHI_6,1))/2
if STATION_BARCODE_EDGES_6[3] >= 180.:  # accounts for wrap around at station 3
    STATION_BARCODE_EDGES_6[3] = STATION_BARCODE_EDGES_6[3]-180.
else:
    STATION_BARCODE_EDGES_6[3] = STATION_BARCODE_EDGES_6[3]+180.
STATION_BARCODE_EDGES_6 = np.append(STATION_BARCODE_EDGES_6, STATION_BARCODE_EDGES_6[0])


# ## Plunging probe calibrations
# Provided by D. Flay.
# Nowhere better to put them.
# Still blinded.

# In[ ]:


# PLUNGING_PROBE_CALIBRATIONS = [-31.574,  # 0
#                                -37.545,  # 1
#                                -27.665,  # 2
#                                -33.719,  # 3
#                                -29.534,  # 4
#                                -13.933,  # 5
#                                -4.149,  # 6
#                                -42.405,  # 7
#                                -40.163,  # 8
#                                -100.103,  # 9
#                                54.498,  # 10
#                                -9.157,  # 11
#                                2.88,  # 12
#                                -40.099,  # 13
#                                -43.37,  # 14
#                                -100.348,  # 15
#                                51.747]  # 16

## pre 28 Sept 2020 calibration
PLUNGING_PROBE_CALIBRATIONS_OLD = [188.537,
                                   182.223,
                                   192.927,
                                   186.147,
                                   190.668,
                                   202.539,
                                   216.653,
                                   175.946,
                                   181.767,
                                   122.287,
                                   274.241,
                                   209.006,
                                   222.794,
                                   175.878,
                                   173.065,
                                   118.348,
                                   268.053]

### unblinded calibration numbers

PLUNGING_PROBE_CALIBRATIONS = [87.521,
                               81.002,
                               91.643,
                               83.039,
                               89.962,
                               103.194,
                               113.698,
                               73.605,
                               80.592,
                               20.797,
                               174.332,
                               107.950,
                               120.301,
                               74.745,
                               72.838,
                               18.166,
                               169.407]
