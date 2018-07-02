# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 15:01:38 2018

@author: Jordy
"""
import sys
sys.path.append('../')
import pandas as pd
import numpy as np
import scipy.interpolate as interpolate
from Aerodynamics import Wing as Awing
from Geometry import Geometry
import math as m


data = pd.read_csv('aerodynamic_data_ms15.dat', ' ', header=None).values

def lookup_data(alpha, ca_c, da):
    # Looksup and interpolates Cl and Cd based on alpha, ca_c and da
    alpha = round(alpha,5)
    alpha = m.degrees(alpha)
    indexca_c = int(100*ca_c)-1
    if da == 0.:
        indexda = abs(da)//5
        localdata = data[int(indexda*50*51+indexca_c*51):int(indexda*50*51+(indexca_c+1)*51),:]
        non_zero_max = max(np.argwhere(localdata[:, 0]))[0]  # last non-zero row
        localdata = localdata[:non_zero_max+1,:]
    else:
        index1da = abs(da)//5
        index2da = abs(da)//5 + 1
        index3da = 0
        localdata1 = data[int(index1da*50*51+indexca_c*51):int(index1da*50*51+(indexca_c+1)*51),:]
        localdata2 = data[int(index2da*50*51+indexca_c*51):int(index2da*50*51+(indexca_c+1)*51),:]
        localdata3 = data[int(index3da*50*51+indexca_c*51):int(index3da*50*51+(indexca_c+1)*51),:]
        non_zero_max1 = max(np.argwhere(localdata1[:, 0]))[0]  # last non-zero row
        non_zero_max2 = max(np.argwhere(localdata2[:, 0]))[0]  # last non-zero row
        non_zero_max3 = max(np.argwhere(localdata3[:, 0]))[0]  # last non-zero row
        non_zero_max = min((non_zero_max1,non_zero_max2,non_zero_max3))
        
        localdata1 = localdata1[:non_zero_max+1,:]
        localdata2 = localdata2[:non_zero_max+1,:]    
        localdata3 = localdata3[:non_zero_max+1,:]   
        
        localdata_alpha = localdata2[:,0]
        localdata_cl    = (localdata2[:,1] - localdata1[:,1])/5 *abs(da) + localdata3[:,1]
        localdata_cd    = (localdata2[:,2] - localdata1[:,2])/5 *abs(da) + localdata3[:,2]
        localdata_cm    = (localdata2[:,4] - localdata1[:,4])/5 *abs(da) + localdata3[:,4]
    
        localdata = np.array([localdata_alpha,
                              localdata_cl,
                              localdata_cd,
                              localdata_cd,
                              localdata_cm])
        localdata = np.transpose(localdata)
    
    Cl_local = interpolate.interp1d(localdata[:,0], localdata[:,1], 'linear', fill_value='extrapolate')
    Cd_local = interpolate.interp1d(localdata[:,0], localdata[:,2], 'linear', fill_value='extrapolate')
    Cm_local = interpolate.interp1d(localdata[:,0], localdata[:,4], 'linear', fill_value='extrapolate')
    if da >= 0:
        Cl = Cl_local(alpha)
        Cd = Cd_local(alpha)
        Cm = Cm_local(alpha)
    else:
        Cl = -Cl_local(-alpha)
        Cd = Cd_local(-alpha)
        Cm = -Cm_local(-alpha)
    
    Cn_lookup = m.cos(m.radians(alpha))*Cl + m.sin(m.radians(alpha))*Cd
    if Cn_lookup !=0:
        xcp = 0.25-Cm/Cn_lookup
    else:
        xcp = 0.25
    return Cl, Cd, Cm, xcp

e = Awing.Oswald_e
AR = Geometry.Wing.A


#Induced angle of attack
alpha_nose = m.radians(5)
alpha_i = 0.
for i in range(10):
    alpha_e = alpha_nose - alpha_i
    Cl, dummy, dummy, dummy = lookup_data(alpha_e,0.1,0)
    alpha_i = Cl / (m.pi * AR * e)
    print(m.degrees(alpha_e))
alpha_nose -= alpha_i
print (alpha_i)

#Induced drag
Cdi = Cl**2 / (m.pi * AR * e)

#Downwash
downwash = 2 * Cl / (m.pi * AR * e)