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
from Geometry import Geometry
from Aerodynamics import Wing as Aero_wing
from Aerodynamics import HT as Aero_HT
from Aerodynamics import VT as Aero_VT
from Aerodynamic import General as Aero_gen
from Performance import Performance


# This file will calculate the ranges for l_h,
# Xlemac, S_h and S_v in which StefX will havel level 1 flying qualities

# Input parameters
gamma_0 = 0                 # Assuming level flight
A = Geometry.Wing.A
CL_alpha = Aero_wing.CL_alpha
MTOW = Geometry.Masses.W_MTOW
S_wing = Geometry.Wing.S
B = Geometry.Wing.S
V_a = Performance.V_a_clean
rho = Performance.rho_a
g0 = Performance.g0
Oswald_e = Aero_wing.Oswald_e
CNH_alpha = Aero_HT.C_Nh_alpha
dE_dalpha = Aero_wing.de_da
Vh_V = Aero_HT.Vh_v
Cbar = Geometry.Wing.MAC
I_yy = 
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