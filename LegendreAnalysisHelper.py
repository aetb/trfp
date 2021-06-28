import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz


import gm2
import trfp
import helper_function_candidates as helper
from numpy import sqrt

## Legendre polynomial moments
legendre_coeffs = np.array([[1./2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # const
                            [0, sqrt(3)/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # x
                            [0, 0, sqrt(3)/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # y
                            [sqrt(5)/-4., 0, 0, sqrt(5)*3./4., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # x^2
                            [0, 0, 0, 0, 3./2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  #xy
                            [sqrt(5)/-4., 0, 0, 0, 0, sqrt(5)*3./4., 0, 0, 0, 0, 0, 0, 0, 0, 0],  # y^2
                            [0, sqrt(7)*-3./4., 0, 0, 0, 0, sqrt(7)*5./4., 0, 0, 0, 0, 0, 0, 0, 0],  # x^3
                            [0, 0, sqrt(15)*-1./4., 0, 0, 0, 0, sqrt(15)*3./4., 0, 0, 0, 0, 0, 0, 0],  #x^2 y
                            [0, sqrt(15)*-1./4., 0, 0, 0, 0, 0, 0, sqrt(15)*3./4., 0, 0, 0, 0, 0, 0], # x y^2
                            [0, 0, sqrt(7)*-3./4., 0, 0, 0, 0, 0, 0, sqrt(7)*5./4., 0, 0, 0, 0, 0],  # y^3
                            [9./16., 0, 0, -90./16., 0, 0, 0, 0, 0, 0, 105./16., 0, 0, 0, 0],  # x^4
                            [0, 0, 0, 0, sqrt(21)*-3./4., 0, 0, 0, 0, 0, 0, sqrt(21)*5./4., 0, 0, 0],  # x^3 y
                            [5./8., 0, 0, -15./8., 0, -15./8., 0, 0, 0, 0, 0, 0, 45./8., 0, 0],  # x^2 y^2
                            [0, 0, 0, 0, sqrt(21)*-3./4., 0, 0, 0, 0, 0, 0, 0, 0, sqrt(21)*5./4., 0],  # x y^3
                            [9./16., 0, 0, 0, 0, -90./16., 0, 0, 0, 0, 0, 0, 0, 0, 105./16.]]  # y^4
                          )

## Cylindrical Legendre polynomial coefficients
## Has a weighting term (711.2/4.5 + x) in dot product
## Not precise due to stupid 711.2/4.5 term, could probably precise-ify it

cyl_leg_coeffs = np.array([[0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # const
                           [-0.0001453, 0.06889, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # x
                           [0, 0, 0.06889, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # y
                           [-0.04447, -0.0002251, 0, 0.1334, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # x^2
                           [0, 0, -0.0002517, 0, 0.1193, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # xy
                           [-0.04447, 0, 0, 0, 0, 0.1334, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # y^2
                           [0.0001427, 0.1578, 0, -0.000428, 0, 0, 0.2631, 0, 0, 0, 0, 0, 0, 0, 0],  # x^3
                           [0, 0, -0.07702, 0, -0.0003899, 0, 0, 0.2311, 0, 0, 0, 0, 0, 0, 0],  # x^2 y
                           [0.0001624, 0.07702, 0, 0, 0, 0.0004873, 0, 0, 0.2311, 0, 0, 0, 0, 0, 0],  # x y^2
                           [0, 0, -0.1578, 0, 0, 0, 0, 0, 0, 0.2631, 0, 0, 0, 0, 0],  # y^3
                           [0.04474, 0.0005033, 0, -0.4474, 0, 0, -0.0008389, 0, 0, 0, 0.522, 0, 0, 0, 0],  # x^4
                           [0, 0,  0.0002471, 0, -0.2734, 0, 0, 0.0007414, 0, 0, 0, 0.4557, 0, 0, 0],  # x^3 y
                           [0.04972, 0.0002517, 0, -0.1491, -0.1491, 0, 0, 0, -0.000755, 0, 0, 0, 0.4474, 0, 0],  # x^2 y^2
                           [0, 0, 0.0005766, 0, -0.2734, 0, 0, 0, 0, -0.000961, 0, 0, 0, 0.4557, 0],  # x y^3
                           [0.04474, 0, 0, 0, 0, -0.4474, 0, 0, 0, 0, 0, 0, 0, 0, 0.522]]  # y^4
                         )
def legendre_moment(l, x, y, cyl=True):
    num = len(x)
    terms = np.array([np.ones(num), x, y, x**2, x*y, y**2, x**3, x**2*y, x*y**2,
                      y**3, x**4, x**3*y, x**2*y**2, x*y**3, y**4])
   
    values = np.empty([num,15])
    for ii in range(15): values[:,ii] = terms[ii]
       
    if cyl: coeffs = cyl_leg_coeffs
    else: coeffs = legendre_coeffs
       
    output = np.matmul(values, coeffs[l])
    output[np.abs(output)<1.e-12] = 0
   
    return output

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

def ThetaL(Multipole, geom, badProbes, tpMomentCap, fpMomentCap, cylindrical=True): 

    lenOrigProbes = len(geom.probes) + len(badProbes)
    
    if (lenOrigProbes > 6):
        momentCap = tpMomentCap
    else:
        momentCap = fpMomentCap
    
    if (len(geom.probes) < momentCap):
        print("Moment cap is higher for the number of functional probes")
    
    _MULTS_FP_N = np.array([Multipole(i, geom.xpos, geom.ypos, cyl=cylindrical)\
                                for i in range(momentCap)])
        
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
        for i in range(lenOrigProbes - momentCap):
            inv_mult = np.concatenate((inv_mult, np.zeros((1,lenOrigProbes))))
        return(inv_mult)
    else:
        for i in range(lenOrigProbes - momentCap):
            INV_MULT = np.concatenate((INV_MULT, np.zeros((1,lenOrigProbes))))
            
        return(INV_MULT)
    
def ProbeDropsL(geom, tpMomentCap = 14, fpMomentCap = 5, cylindrical = True):
    reducedGeometry, badProbes = DropPos(geom)
    
    return ThetaL(legendre_moment, reducedGeometry, badProbes, tpMomentCap, fpMomentCap, cylindrical=cylindrical)

#############################################################
###############Start of 3D Legendre stuff####################
#############################################################

#Make a class similar to the Geometry file, made for dropping probes/ drop with 2 options:geometrical or number
def legendre_moment3d(LNS, x, y, phi, cyl = True):
    """ Computes 3D Legendre Moments.
    
    Takes regular Legendre Moment as developed by Alec Tewsley-Booth and 
    adds a cos(n \phi) or sin(n\phi term). L corresponds to the Legendre
    Moment order (L0,...,Lm), N corresponds to the azimuthal order or
    periodicity (0 ,..., n). Cosine is "normal" while sine is "skew". Note
    that are is no zeroth order skew Legendre Moments since Sin(0) = 0.
    
    Args:
    LNS: a (3, X) numpy array containing the moments to be calculated;
         X is any length > 1 this array is usually provided by the
         FlattenMoment function.
    
    x:   a 1-d array containing the x-positions of the P points where the
         Legendre Moments are to be calculated.
         
    y:   a 1-d array containing the y-positions of the P points where the
         Legendre Moments are to be calculated.
         
    phi: a 1-d array containing the phi-positions (in radians) of the P
         points where the Legendre Moments are to be calculated.
         
    cyl: argument that if True calculates cylindrical Legendre Moments. If
         False, calculates cartesian Legendre Moments.
    
    """
    s = LNS[2]
    if s == 1:
        return legendre_moment(LNS[0],x,y,cyl)*np.sin(LNS[1]*phi)
    elif s == 0:
        return legendre_moment(LNS[0],x,y,cyl)*np.cos(LNS[1]*phi)
    
class RingGeometry:
    """
    A class that stores the x, y and phi positions of the probes to be used.
    In the case of fixed probes, it also containes which probes are function-
    al. Additional functions to simplify other functions are also here. x, y,
    phi and p have the same format as trfp.STATION_PROBE_ID for fixed probes.
    """
    def __init__(self, x, y, phi, p):
        """
        Initialized the class, not much to say here except that p is an array
        of 0s and 1s where 0 means that such a probe is non-functional and
        1 means it is functional
        """
        self.xpos = x
        self.ypos = y
        self.phis = phi
        self.probes = p
    
    def IsReduced(self):
        """
        Checks to see whether a geometry instance is reduced. If a geometry
        class is reduced then all probes are functional. Otherwise it is not
        reduced.
        
        These are some tests:
        >>> x = [np.arange(5)]
        >>> y = [np.arange(5)]
        >>> phi = [np.linspace(0, 360, 5)]
        >>> p = [np.zeros(5) +1]
        >>> geom = RingGeometry(x ,y ,phi ,p)
        >>> geom.IsReduced()
        True
        
        >>> p = [np.array([0,1,1,1,1])]
        >>> geom = RingGeometry(x ,y ,phi ,p)
        >>> geom.IsReduced()
        False
        """
        isReduced = True
        
        for i in range(len(self.probes)):
            for j in range(len(self.probes[i])):
                if self.probes[i][j] != 1:
                    return False
        return isReduced
    
    def getBadProbes(self):
        badProbes = []
        for i in range(len(self.probes)):
            badProbesInThisStation = np.array([])
            for j in range(len(self.probes[i])):
                if not self.probes[i][j]:
                    badProbesInThisStation = np.append(badProbesInThisStation,j)
                
            badProbes.append(badProbesInThisStation.astype(int))
        return badProbes
    
    def getNumberOfBadProbes(self):
        '''
        
        '''
        N = 0
        badProbes = self.getBadProbes()
        for i in range(len(badProbes)):
            for j in range(len(badProbes[i])):
                N+=1
        return N
                
    
    def getNumberOfProbes(self):
        N = 0
        for i in range(len(self.probes)):
            N+= len(self.probes[i])
        return N

def Drop3dPos(geom): #can make XYPhiProbeArray
    badProbes = geom.getBadProbes()
    newx = []
    newy = []
    newPhis = []
    newProbes = []
    
    for i in range(len(geom.probes)):
        if len(badProbes[i]) == len(geom.probes[i]):
            continue
        PROBES = geom.probes[i].copy()
        X = geom.xpos[i].copy()
        Y = geom.ypos[i].copy()
        PHI = geom.phis[i].copy()
        
        newx.append(np.delete(X, badProbes[i]))
        newy.append(np.delete(Y, badProbes[i]))
        newPhis.append(np.delete(PHI, badProbes[i]))
        newProbes.append(np.delete(PROBES, badProbes[i]))
    ReducedGeometry = RingGeometry(newx, newy, newPhis, newProbes)
    return ReducedGeometry, badProbes

def FlattenMoments(lcap, ncap):
    #First lcap Lmoments and first ncap Az. moments (n=0 counts as 1 too)
    Ls = np.tile(np.repeat(np.arange(lcap),2), ncap)
    Ns = np.sort(np.repeat(np.arange(ncap), lcap*2))
    Ss = np.tile([0,1], lcap*ncap)
    return np.delete(Ls,np.arange(lcap)*2 +1),\
           np.delete(Ns,np.arange(lcap)*2 +1),\
           np.delete(Ss,np.arange(lcap)*2 +1)

def CreateLNSHeader(lcap,ncap): #n1s,l3
    Ls, Ns, Ss = FlattenMoments(lcap, ncap)
    Header = []
    for i in range(len(Ls)):
        col = 'N' + str(Ns[i])
        if Ss[i] == 1:
            col = col + 's,'
        else:
            col = col + 'n,'
        col = col + 'l' + str(Ls[i])
        Header.append(col)
    return Header

def InvMultipole(MomentFxn, geom, lcap, ncap, cylindrical = True):
    Xs = np.concatenate(geom.xpos)
    Ys = np.concatenate(geom.ypos)
    Phis = np.concatenate(geom.phis)
    
    LNS = np.zeros((2*lcap*ncap - lcap, 3), dtype = int)
    LNS[:,0], LNS[:,1], LNS[:,2] = FlattenMoments(lcap, ncap)
    
    
    if (len(LNS)>geom.getNumberOfProbes()):
        print("Not enough probes for the Moments, overfitting...")
        
    MULTS = np.array([MomentFxn(LNS[i], Xs, Ys, Phis, cyl = cylindrical)\
                                for i in range(len(LNS))])
    
    MULTS[np.abs(MULTS) < 1.0e-9] = 0
    
    INV_MULT = np.linalg.pinv(np.transpose(MULTS))
    INV_MULT[np.abs(INV_MULT) < 1.0e-9] = 0
    
    I = np.dot(INV_MULT, np.transpose(MULTS))
    I[np.abs(I) < 1.0e-9] = 0
    print(np.sum(np.abs(I)) - min(len(I), len(I[0])))
    
    return INV_MULT
    
    
#put this above ProbeDropsL Also, must make a copy of badProbes
def Theta3dL(MomentFxn, geom, badProbes, lcap, ncap, cylindrical = True):
    lenOrigProbes = 378
    INV_MULT = InvMultipole(MomentFxn, geom, lcap, ncap, cylindrical)
    
    inv_mult = np.zeros((len(INV_MULT), lenOrigProbes))
    k = 0
    l = 0 #This counts how many iterations have been done in the loop
    if geom.getNumberOfProbes() < lenOrigProbes:
        for i in range(len(trfp.STATION_PROBE_ID)):
            for j in range(len(trfp.STATION_PROBE_ID[i])):
                if (len(badProbes[i]) > 0) and (j == badProbes[i][0]):
                    badProbes[i] = np.delete(badProbes[i], [0]) #"skip"
                    l+=1
                else:
                    inv_mult[:,l] = INV_MULT[:,k]
                    k+=1
                    l+=1
        return inv_mult
    else:
        return INV_MULT

def ProbeDrops3dL(MomentFxn, geom, lcap, ncap, cylindrical = True):
    reducedGeometry, badProbes = Drop3dPos(geom)
    
    return Theta3dL(MomentFxn, reducedGeometry, badProbes, lcap, ncap, cylindrical=cylindrical)

######################################################################
########################End of 3D Leg Stuff###########################
######################################################################


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
                   fixedProbes = np.zeros((76,6))+1, tpMomentCap = 14,\
                   fpMomentCaps = np.zeros(len(trfp.STATION_PROBE_ID),dtype='int')+4):
    #Make regular regular 4 and 6 FP thetas -JS
    FP6Geometry = Geometry(trfp.FP6_X, trfp.FP6_Y, np.zeros(6)+1)
    FP6Theta = ProbeDropsL(FP6Geometry, 14, 5) #14 trm and m5 for fp
        
    FP4Geometry = Geometry(trfp.FP4_X, trfp.FP4_Y, np.zeros(4)+1)
    FP4Theta = ProbeDropsL(FP4Geometry, 14, 4) #same as FP6Theta
    
    tr_run = 'tr_phi' in interp_df.columns
    
    moment_df = pd.DataFrame(index=interp_df.index)
    
    # calculate trolley moments if needed
    if tr_run:
        #Trolley Probe Dropping -JS
        TrolleyGeometry = Geometry(trfp.TR_X, trfp.TR_Y, trolleyProbes)
        TrolleyTheta = ProbeDropsL(TrolleyGeometry, tpMomentCap=tpMomentCap)
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
        if (len(trfp.STATION_PROBE_ID[station]) == 6):
            if not(0 in fixedProbes[station]) and (fpMomentCaps[station] == 5):
                theta_fp = FP6Theta
                #theta_fp = trfp.THETA_FP_6
            else: 
                StGeom = Geometry(trfp.FP6_X, trfp.FP6_Y, fixedProbes[station])
                StTheta = ProbeDropsL(StGeom, fpMomentCap=fpMomentCaps[station])
                theta_fp = StTheta
        
        elif(station == 41):
            St41Geom = Geometry(trfp.FP4_X_ST41, trfp.FP4_Y, fixedProbes[station][:4])
            St41Theta = ProbeDropsL(St41Geom, fpMomentCap=fpMomentCaps[station])
            theta_fp = St41Theta
        
        elif (station == 37) or (station == 39):
            St37_39Geom = Geometry(trfp.FP4_X_ST37_ST39, trfp.FP4_Y,\
                                   fixedProbes[station][:4])
            St37_39Theta = ProbeDropsL(St37_39Geom, fpMomentCap=fpMomentCaps[station])
            theta_fp = St37_39Theta
            
        else:
            if (not (0 in fixedProbes[station][:4])) and (fpMomentCaps[station] == 5):
                theta_fp = FP4Theta
            else:
                StGeom = Geometry(trfp.FP4_X, trfp.FP4_Y, fixedProbes[station][:4])
                StTheta = ProbeDropsL(StGeom, fpMomentCap=fpMomentCaps[station])
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