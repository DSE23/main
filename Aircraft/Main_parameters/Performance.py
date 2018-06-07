"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

StefX Flight performance parameters
"""

import sys
import math as m
import numpy as np
sys.path.append('../')

from Control import Calc_ISA_100km as ISA

mtow = 889.8                    # [kg] Max. Take-off weight
cl_max_hld = 2.522              # [-] CL max with HLD's deployed
cl_max_clean = 1.3154           # [-] CL max in clean config
to_distance = 400              # [m] Take-off distance
top = 130                       # [-] Take-Off parameter
s_land = 550                    # [m] Landing distance
V_cruise = 94.1                 # [m/s] Cruise speed
h_c = 3000                      # [m] Cruise altitude !!!please check!!!
rho_c = ISA.isacal(h_c)[2]      # Cruise density
Temp_c = ISA.isacal(h_c)[1]     # Cruise temp

