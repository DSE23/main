"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

StefX Flight performance parameters
"""

import sys
import math as m
import numpy as np
sys.path.append('../')
from Geometry import Geometry
from Misc import Q_, ureg
from Control import Calc_ISA_100km as ISA
from Aerodynamics import Wing as Aero_wing

mtow = Geometry.Masses.W_MTOW       # [kg] Max. Take-off weight
cl_max_hld = Aero_wing.CL_max_hld   # [-] CL max with HLD's deployed
cl_max_clean = Aero_wing.CL_max_hld # [-] CL max in clean config
to_distance = Q_("400 m")              # [m] Take-off distance
top = 130                      # [-] Take-Off parameter
s_land = Q_("550 m")                   # [m] Landing distance
V_cruise = Q_("94.1 m/s")                 # [m/s] Cruise speed
h_c = Q_("3000 m")                     # [m] Cruise altitude !!!please check!!!
rho_c = ISA.isacal(h_c.magnitude)[2]   # Cruise density
rho_c *= ureg("kg/(m**3)")
Temp_c = ISA.isacal(h_c.magnitude)[1]     # Cruise temp
Temp_c *= ureg("K")
gamma_isa = 1.41                # Specific heat ratio
m_c = (V_cruise/np.sqrt(gamma_isa*Temp_c*rho_c)).magnitude  # [-] Cruise Mach number
roc = Q_("16 m/s")                     # [m/s] Rate of Climb
Climb_angle = Q_("45 deg")                # [deg] Climb angle
V_sust = Q_("67.135 m/s")                # [m/s] Sustained turn velocity
n_sust = 4.16                   # [-] Sustained turn load factor
h_a = Q_("600 m")                      # [m] Manouevring height
rho_a = ISA.isacal(h_a.magnitude)[2]       # [kg/m^3] Manouevring density
rho_a *= Q_(" 1 kg/m**3")
S_wing = Geometry.Wing.S                 # Wing span for next calc.
g0 = Q_("9.80665 m/s**2")
V_a_hld = np.sqrt((2*g0*mtow)/(cl_max_hld*rho_a*S_wing))   # Man. Veloc. at with HLD
V_a_clean = np.sqrt((2*g0*mtow)/(cl_max_clean*rho_a*S_wing))   # Man. Veloc. clean

