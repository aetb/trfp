#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd
import scipy

from scipy.optimize import curve_fit


# In[4]:


input_12 = np.loadtxt('/home/aetb/pythonScripts/trfp/mu_avg/sample_st12.txt')
input_18 = np.loadtxt('/home/aetb/pythonScripts/trfp/mu_avg/sample_st18.txt')

x_12 = input_12[:,0].reshape(60,60).transpose()
y_12 = input_12[:,1].reshape(60,60).transpose()
N_12 = input_12[:,2].reshape(60,60).transpose()

x_18 = input_18[:,0].reshape(60,60).transpose()
y_18 = input_18[:,1].reshape(60,60).transpose()
N_18 = input_18[:,2].reshape(60,60).transpose()

x = np.arange(-59,60,2)
y = np.arange(-59,60,2)

x_dist_12 = np.sum(N_12,axis=0)
x_dist_12 /= np.max(x_dist_12)
y_dist_12 = np.sum(N_12,axis=1)
y_dist_12 /= np.max(y_dist_12)

x_dist_18 = np.sum(N_18,axis=0)
x_dist_18 /= np.max(x_dist_18)
y_dist_18 = np.sum(N_18,axis=1)
y_dist_18 /= np.max(y_dist_18)

def _x_dist_long(x, A1, A2, A3, x1, x2, x3, w1, w2, w3):
    f1 = A1 * np.exp(-(x-x1)**2/2/w1**2)
    f2 = A2 * np.exp(-(x-x2)**2/2/w2**2)
    f3 = A3 * np.exp(-(x-x3)**2/2/w3**2)
    return f1+f2+f3

def _y_dist(y, A, y0, wy):
    f1 = np.exp(-(y-y0)**2/2/wy**2)
    return A * f1 / wy / np.sqrt(2*np.pi)

coeffs_x_12, _ = curve_fit(_x_dist_long, x, x_dist_12, p0=[0.17, 0.5, 0.9, -24, -2, 20, 8.3, 12.7, 9.9])
coeffs_x_18, _ = curve_fit(_x_dist_long, x, x_dist_18, p0=[0.17, 0.5, 0.9, -23, -2, 20, 8.3, 12.7, 9.9])

coeffs_y_12, _ = curve_fit(_y_dist, y, y_dist_12, p0=[1, 30, 10])
coeffs_y_18, _ = curve_fit(_y_dist, y, y_dist_18, p0=[1, 30, 10])

## need to find averages, even if they are out of order
## order by center paramters
x_order_12 = np.argsort(coeffs_x_12[3:6])
x_order_18 = np.argsort(coeffs_x_18[3:6])

x23 = ((coeffs_x_12[3:6][x_order_12][1] - coeffs_x_12[3:6][x_order_12][2]
        + coeffs_x_18[3:6][x_order_18][1] - coeffs_x_18[3:6][x_order_18][2])
       / (coeffs_x_12[3:6][x_order_12][0] - coeffs_x_12[3:6][x_order_12][2]
          + coeffs_x_18[3:6][x_order_18][0] - coeffs_x_18[3:6][x_order_18][2]))

w1 = np.abs(coeffs_x_12[6:9][x_order_12][0]/(coeffs_x_12[3:6][x_order_12][0] - coeffs_x_12[3:6][x_order_12][2])
      + coeffs_x_18[6:9][x_order_18][0]/(coeffs_x_18[3:6][x_order_18][0] - coeffs_x_18[3:6][x_order_18][2])
     )/2
w2 = np.abs(coeffs_x_12[6:9][x_order_12][1]/(coeffs_x_12[3:6][x_order_12][0] - coeffs_x_12[3:6][x_order_12][2])
      + coeffs_x_18[6:9][x_order_18][1]/(coeffs_x_18[3:6][x_order_18][0] - coeffs_x_18[3:6][x_order_18][2])
     )/2
w3 = np.abs(coeffs_x_12[6:9][x_order_12][2]/(coeffs_x_12[3:6][x_order_12][0] - coeffs_x_12[3:6][x_order_12][2])
      + coeffs_x_18[6:9][x_order_18][2]/(coeffs_x_18[3:6][x_order_18][0] - coeffs_x_18[3:6][x_order_18][2])
     )/2

A13 = (coeffs_x_12[0:3][x_order_12][0]/coeffs_x_12[0:3][x_order_12][2]
       + coeffs_x_18[0:3][x_order_18][0]/coeffs_x_18[0:3][x_order_18][2]
      )/2
A23 = (coeffs_x_12[0:3][x_order_12][1]/coeffs_x_12[0:3][x_order_12][2]
       + coeffs_x_18[0:3][x_order_18][1]/coeffs_x_18[0:3][x_order_18][2]
      )/2

def _x_dist(x, A, x3, x31):
    f3 = np.exp(-(x-x3)**2/2/(w3*x31)**2)
    f1 = A13 * np.exp(-(x-x3+x31)**2/2/(w1*x31)**2)
    f2 = A23 * np.exp(-(x-x3+x23*x31)**2/2/(w2*x31)**2)
    
    F3 = np.sqrt(2*np.pi) * w3 * x31
    F1 = A13 * np.sqrt(2*np.pi) * w1 * x31
    F2 = A23 * np.sqrt(2*np.pi) * w2 * x31
    
    return A*(f1 + f2 + f3)/(F1+F2+F3)


# In[5]:


def _xy_dist_fit(XY, A, x3, x31, y0, wy):
    x, y = XY
    fx = _x_dist(x, 1, x3, x31)
    fy = _y_dist(y, 1, y0, wy)
    return A * fx * fy

def _xy_dist(x, y, A, x3, x31, y0, wy):
    fx = _x_dist(x, 1, x3, x31)
    fy = _y_dist(y, 1, y0, wy)
    return A * fx * fy


# In[6]:


## include azimuthal dependence

R0 = 7112

def _x3(theta):  # just a sample, would need to be filled by radial closed orbit
    return 19.4 + 0.05*np.cos(theta)
def _x31(theta):  # needs to be fit by betatron breathing?
    return 45.9 * (1 + 0.05*np.cos(4*theta))
def _y0(theta):  # vertical offset/closed orbit
    return -0.17
def _wy(theta):
    return 13.5 * (1 - 0.1*np.cos(4*theta))

def _xyz_dist(x, y, theta, A):
    x3 = _x3(theta)
    x31 = _x31(theta)
    y0 = _y0(theta)
    wy = _wy(theta)
    return _xy_dist(x, y, 1, x3, x31, y0, wy) / (2*np.pi*R0)  # crappy normalization

def _xyz_integrand(x, y, theta, A):
    return (R0+x) * _xyz_dist(x, y, theta, A)


# In[ ]:




