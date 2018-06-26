""" This is a tool to calculate for every speed the required power and AoA to maintain a certain climb angle"""

import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
import time
from Misc import ureg, Q_
from Geometry import Geometry as GM
from Aerodynamics import Aeroprops as Aeroprops
from Aerodynamics import Wing as AWing
import Performance as PF
from Misc import Init_parm as IP
import matplotlib.pyplot as plt
import math as m
from Propulsion_and_systems import Propdata as Propdata

# Import parameters
P_to = PF.P_to
C_d_0 = Aeroprops.CD0_tot
mass = GM.Masses.W_MTOW
W = mass * Q_("9.81 m/s**2")
S = GM.Wing.S
rho = Q_("1.225 kg/(m**3)")
C_L_alpha = AWing.CL_alpha
C_L_alpha *= Q_("1/rad")
C_l_max = AWing.CL_max
A = GM.Wing.A
e = AWing.Oswald_e
eta_prop = PF.eta_prop
dp = PF.dp
V_stall =  PF.V_stall_hld
alpha_max = AWing.alpha_stall
V_lof = 1.05 * V_stall
V_max = PF.V_a_clean
#Tmax = (P_to**2*eta_prop**2*m.pi*dp**2/2*rho)**(1/3)


# Inputs
flight_path_angle = Q_("20 deg")
dV = Q_("1 m/s")


# Empty lists
V_list =[]
alpha_list = []
T_list = []
T_max_list =[]


flight_path_angle.ito(Q_("rad"))


for V in np.arange(16, 100, 1):
    V *= Q_("m/s")

    L_req = W * np.cos(flight_path_angle)
    alpha_req = L_req / (0.5 * rho * V**2 * S * C_L_alpha)
    C_D = C_d_0 + (C_L_alpha * alpha_req)**2/ (m.pi * A * e)
    D = C_D * 0.5 * rho * V**2 *S
    T_req = D + W * np.sin(flight_path_angle)
    alpha_req.ito(Q_("deg"))

    V_list.append(V.magnitude)
    alpha_list.append(alpha_req.magnitude)
    T_list.append(T_req.magnitude)
    T_max = Propdata.data_reader(V, "Total thrust")
    T_max_list.append(T_max.magnitude)

    #print(V, alpha_req, P_ratio)


P_to.ito(Q_("hp"))
plt.figure(1)
plt.plot(V_list, T_list)
plt.plot(V_list, T_max_list, color ='r')
plt.xlabel('Velocity  [m/s]')
plt.ylabel('Thrust needed [N]')
#plt.axhline(y= P_to.magnitude, color='r')

plt.figure(2)
plt.plot(V_list, alpha_list, )
plt.xlabel('Velocity  [m/s]')
plt.ylabel('AoA needed [rad]')
plt.axhline(y= alpha_max, color='r')

plt.show()
