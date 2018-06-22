"""
Name: Fuel Tanks
Department: Propulsion and Aircraft Systems
Last updated: 22/06/2018 14:17 by Ties
"""

import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

import numpy as np
from Geometry import Geometry
from Aerodynamics import Aeroprops
from Performance import Performance
from Propulsion_and_systems import Propdata

# Find CL & CD for max Cl/Cd (cruise)

C_L = np.sqrt(np.pi*Geometry.Wing.A*Aeroprops.CD0_tot)
C_D = 2*Aeroprops.CD0_tot

# Define weights
W_fuel = Q_("100 kg")
W_empty = Geometry.Masses.W_OEW + Geometry.Masses.W_pilot + W_fuel

# Calculate Drag
D = C_D/C_L * (W_empty+W_fuel)*Performance.g0
D.ito(ureg.newton)
print("Drag: {}".format(D))

# Calculate V_cruise
V_cruise = np.sqrt(((W_empty + W_fuel)*Performance.g0) / (0.5*Performance.rho_c*C_L*Geometry.Wing.S))
print("Cruise velocity: {}".format(V_cruise))

# Get power needed from prop file
P = Propdata.data_reader(V_cruise, "Power needed")
P.ito(ureg.hp)
print("Power needed: {}".format(P))