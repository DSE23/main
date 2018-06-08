# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

Stability derivatives range for level 1 flying qualities
"""

import sys
import math as m
import numpy as np
sys.path.append('../')

from Misc import Init_parm as IP
from Geometry import Wing

# This file will calculate the ranges for l_h,
# Xlemac, S_h and S_v in which StefX will havel level 1 flying qualities

# Input parameters
gamma_0 = 0                 # Assuming level flight
A = Wing.A
CL_alpha = IP.CL_alpha
MTOW = IP.MTOW
S_wing = IP.S_wing
B = Wing.s
V_a = IP.V_a_clean
rho = IP.rho0
g0 = IP.g0
Oswald_e = IP.Oswald_e
CNH_alpha = IP.CN_h_alpha
dE_dalpha = IP.Downwash
Vh_V = IP.Vh_V
Cbar = IP.MAC
I_yy = IP.I_yy
g = IP.g0
CNW_alpha = IP.CN_w_alpha
# specific parameters (only check if either Lambda or A changes)


# Iteration values
S_h = 2.629
l_h = 5.2
# Calculated values
C_L = MTOW/(0.5*rho*V_a**2*S_wing)


# Stability Derivatives

CX0 = (MTOW * g)/(0.5 * rho * V_a**2 * S_wing) * m.sin(gamma_0)
CXu = -2 * C_L * m.tan(gamma_0)
CX_alpha = C_L * (1-(2*CL_alpha))/(Oswald_e*A*m.pi)
CZ0 = -(MTOW * g)/(0.5*rho*V_a**2*S_wing)*m.cos(gamma_0)
CZu = -2*C_L
CZ_alpha = -CL_alpha
CZ_alphadot = -CNH_alpha * Vh_V**2 * dE_dalpha * S_h * l_h / (S_wing * Cbar)
CZq = -2 * Vh_V**2 * S_h * l_h / (S_wing * Cbar)
Cmu = 0
Cm_alpha =  CNW_alpha * 