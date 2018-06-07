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
S_wing = 11.74               # Wing surface
W = 822                      # MTOW, updates from structure
rho0 = 1.225                    # Density at sea level
g0 = 9.80665
V_stall = m.sqrt((2*W*g0)/(rho0*S_wing*CLmax))
n_max = 10                      # Max load factor
V_man = V_stall*m.sqrt(n_max)

d_delta_a= 60 / 180 *m.pi #rad delta aileron form minus to plus
d_s_a= 0.41 #m      stick deflection
Sa= 2.33401 #m2     aileron surface area
ca= 0.359227 #m     aileron cord
C_h_alpha= 0.009511     #hinge moment
C_l_delta_a= 0.833
C_l_p= -0.632399
C_h_delta= 0.196066
b= 7.541
maxdefl_aileron = 30/180 *m.pi #rad
V_stall = 32 #m/s
V_max= 95 #m/s
alpha_max= 18 #deg
Fa= 86 #Newton stick force
