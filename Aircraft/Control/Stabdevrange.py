# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

Stability derivatives range for level 1 flying qualities
"""

import sys
sys.path.append('../')

from Misc import Init_parm as IP
from Geometry import Wing

# This file will calculate the ranges for l_h,
# Xlemac, S_h and S_v in which StefX will havel level 1 flying qualities

# Input parameters

A = Wing.A
CL_alpha = IP.CL_alpha
MTOW = IP.MTOW
S_wing = IP.S_wing
B = Wing.s
V_a = IP.V_man
rho = IP.rho0
g0 = IP.g0
Oswald_e = IP.Oswald_e
CNH_alpha = 4
# specific parameters (only check if either Lambda or A changes)


#Iteration values
S_h = 2.629

# Calculated values
C_L = MTOW/(0.5*rho*V_a**2*S_wing)


# Stability Derivatives




