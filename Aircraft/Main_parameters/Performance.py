"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

StefX Flight performance parameters
"""

import sys
import math as m
sys.path.append('../')
from Main_parameters import Wing
from Control import Calc_ISA_100km as ISA

mtow = 889.8                    # [kg] Max. Take-off weight
cl_max_hld = 2.522              # [-] CL max with HLD's deployed
cl_max_clean = 1.3154           # [-] CL max in clean config
to_distance = 400               # [m] Take-off distance
top = 130                       # [-] Take-Off parameter
s_land = 550                    # [m] Landing distance
V_cruise = 94.1                 # [m/s] Cruise speed
h_c = 3000                      # [m] Cruise altitude !!!please check!!!
rho_c = ISA.isacal(h_c)[2]      # Cruise density
Temp_c = ISA.isacal(h_c)[1]     # Cruise temp
gamma_isa = 1.41                # You know this gamma
m_c = V_cruise/m.sqrt(gamma_isa*Temp_c*rho_c)  # [-] Cruise Mach number
roc = 16                        # [m/s] Rate of Climb
Climb_angle = 45                # [deg] Climb angle
V_sust = 67.135                 # [m/s] Sustained turn velocity
n_sust = 4.16                   # [-] Sustained turn load factor
h_a = 600                       # [m] Manouevring height
rho_a = ISA.isacal(h_a)[2]      # [kg/m^3] Manouevring density
S_wing = Wing.S                 # Wing span for next calc.
g0 = 9.80665
V_a_hld = m.sqrt((2*g0*mtow)/(cl_max_hld*rho_a*S_wing))   # Man. Veloc. at with HLD
V_a_clean = m.sqrt((2*g0*mtow)/(cl_max_clean*rho_a*S_wing))   # Man. Veloc. clean

