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
#from Aerodynamics import Drag
from Aerodynamics import Wing
from Performance import Performance
from Propulsion_and_systems import Propeller
import matplotlib.pyplot as plt
# Find CL & CD for max Cl/Cd (cruise)

# Calculate Acro fuel needed

fuel_cons_a = Q_("142 lbs/hr")
fuel_density = Q_("0.721 kg/l")
fuel_cons_a = fuel_cons_a/fuel_density
fuel_cons_a.ito("l/hr")
# print(fuel_cons_a)
min_fuel_cap_a = fuel_cons_a/2
print("Minimum amount of acro fuel: {}".format(min_fuel_cap_a))

# C_L = np.sqrt(np.pi*Geometry.Wing.A*Aeroprops.CD0_tot)
# C_D = 2*Aeroprops.CD0_tot
C_D = Aeroprops.CD0_tot

# Define weights
W_fuel = Q_("100 kg")
W_empty = Geometry.Masses.W_OEW + Geometry.Masses.W_pilot + W_fuel
W_tot = (Geometry.Masses.W_OEW + Geometry.Masses.W_pilot + W_fuel) * Performance.g0

# Calculate Drag
# D = C_D/C_L * (W_empty+W_fuel)*Performance.g0
# D.ito(ureg.newton)
# print("Drag: {}".format(D))
#
# # Calculate V_cruise
# V_cruise = np.sqrt(((W_empty + W_fuel)*Performance.g0) / (0.5*Performance.rho_c*C_L*Geometry.Wing.S))
# print("Cruise velocity: {}".format(V_cruise))
#
# # Get power needed from prop file
# P = Propeller.Powercalc(V_cruise.magnitude, D.magnitude)[3]
# P = P*Q_("W")
# P.ito(ureg.hp)
# print("Power needed: {}".format(P))

# Use "economic cruise" engine power
P = Q_("192 hp")
V_cruise = Q_("86.5 m/s")

T1 = Q_("0 newton")
# T2 = []
T2 = Q_("1000 newton")
vlst = []

P.ito(ureg("watt"))
# print(P)

# for i in range(0, 25, 1):
#     V_cruise = 71+i
#     T = Propeller.Thrustcalc(V_cruise)[1]
#     T2.append(T)
#     vlst.append(V_cruise)
#     print(V_cruise, T)
#
# plt.plot(vlst, T2)
# plt.show()

while abs(T1-T2) > Q_("10 newton"):
    if T1 < T2:
        V_cruise += Q_("0.1 m/s")  # abs(T1.magnitude-T2.magnitude)*Q_("m/s")
    elif T1 > T2:
        V_cruise -= Q_("0.1 m/s")  # abs(T1.magnitude-T2.magnitude)*Q_("m/s")

    T1 = (C_D +
          (W_tot / (0.5 * Performance.rho_c * V_cruise**2 * Geometry.Wing.S))**2
           / (np.pi*Wing.Oswald_e*Geometry.Wing.A)) * 0.5*Performance.rho_c*V_cruise**2*Geometry.Wing.S

    eff = 89  # Propeller.ThrustcalcV2(V_cruise.magnitude, P.magnitude)[0]
    T2 = P/V_cruise * (eff/100)
    # T2 = Propeller.Thrustcalc(V_cruise.magnitude)[1]
    # T2 = T2*Q_("newton")
    T2.ito(ureg.newton)
    T1.ito(ureg.newton)

    # print("V: {}".format(V_cruise))
    # print("eff: {}".format(eff))
    # print("T1: {}".format(T1))
    # print("T2: {}".format(T2))

print("V_cruise: {}".format(V_cruise))

fuel_cons_c = Q_("14 gallons / hour")
fuel_cons_c.ito(ureg("liters /hour"))
distance_per_liter = V_cruise/fuel_cons_c
distance_per_liter.ito(ureg("kilometers/liter"))
print("Kilometers per liter: {}".format(distance_per_liter))
min_range = Q_("750 km")
min_fuel_cap = min_range/distance_per_liter
print("Minimum fuel volume needed: {}".format(min_fuel_cap))
