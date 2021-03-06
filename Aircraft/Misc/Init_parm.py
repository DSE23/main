# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:19:31 2018

@author: jurian
"""

import sys
import math as m
import scipy.constants as con
sys.path.append('../')    # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
# This file contains initial parameters found in the Midterm
# If there is a file where a value of this list is calculated
# please link to it, to make sure that all other scripts will
# continue to work, and will use the latest version of the value
# !!!Check Regurlary!!!


CLmax = 2.52                   # CLmax with HLD, should be from aero
CLmax_clean = 1.3154           # From aero guys
CL_alpha = 4.6531              # CL_alpha 1/rad
S_wing = 11.74               # Wing surface
MTOW = 889                      # MTOW, updates from structure
rho0 = 1.225                    # Density at sea level
g0 = 9.80665
V_stall = m.sqrt((2*MTOW*g0)/(rho0*S_wing*CLmax)) #Stall speed not clean condition
V_stall_clean = m.sqrt((2*MTOW*g0)/(rho0*S_wing*CLmax_clean)) #Stall speed clean condition
n_max = 10                      # Max load factor
V_man = V_stall*m.sqrt(n_max)   # Manoeuvre speed @Sea level
V_a_clean = V_stall_clean*m.sqrt(n_max) #Manoeuvre speed clean condition
Oswald_e = 0.79                 # Oswald efficiency factor
CN_h_alpha = 3.238969           # From aero guys
Downwash = 0.7957               # depsilon/dalpha from aero guys
Vh_V = 0.93                     # Vh over V also from aero
I_yy = 1565                     # MMoI from structures
CN_w_alpha = 4.653              # From aero again
MAC = 1.53                      # Mean Aerodynamic chord
S_h = 2.629                     # Horizontal tail surface
X_htail = 5.27                  # Distance 0.25c Htail to nose
X_cg_body_group = 0             # !!!Can somebody calculate this
W_body_group = 612              # kg from structures
W_wing = 121                    # kg from structures
X_cg_wing = 0.5 * MAC           # Distance from Xlemac
W_pilot = 100                   # kg 
W_fuel = 57                     # kg
X_cg_pilot = 2.24               # !!!Please check!!! Distance from nose
X_cg_fuel = 1.18                  # From nose
Xlemac = 1.24                   # Will be iterated for, distance from nose




d_delta_a= 60 / 180 * m.pi #rad delta aileron form minus to plus
d_s_a= 0.41 #m      stick deflection
Sa= 2.33401 #m2     aileron surface area
ca= 0.359227 #m     aileron cord
C_h_alpha= 0.009511     #hinge moment
C_l_delta_a= 0.833
C_l_p= -0.632399
C_h_delta= 0.196066
b= 7.541
maxdefl_aileron = 30/180 *m.pi #rad
alpha_max= 18 #deg
Fa= 86 #Newton stick force
