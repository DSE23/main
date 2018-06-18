""""This file calculates the maximum sustained climb possible at each speed. It takes into account the maximum power available and stall angle.
It assumes that thrust is along the aerodynamic frame"""

import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
from Misc import ureg, Q_
from Geometry import Geometry as GM
from Aerodynamics import Aeroprops as Aeroprops
from Aerodynamics import Wing as AWing
import Performance as PF
from Misc import Init_parm as IP
import matplotlib.pyplot as plt
import math as m



# Get parameters
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
V_stall =  PF.V_stall_clean
alpha_max = AWing.alpha_stall
alpha_max *= Q_("deg")
V_a = PF.V_a_clean
Tmax = (P_to**2*eta_prop**2*m.pi*dp**2/2*rho)**(1/3)

# Inputs
V_step = 1 # m/s
y_step = 0.1 # deg

# make empty lists
V_list = []
y_max_list = []


# The process will be run for all speeds between Vstall and Va
for V in np.arange (V_stall.magnitude ,V_a.magnitude, V_step):
    V *= Q_("m/s")

    # At a new V these lists have to be made empty
    T_D_list =[]
    T_list = []
    D_list = []
    y_list = []

    # For every climb angle the thrust and AoA required are calculated as well as the corresponding drag.
    # If the thrust or AoA is higher than possible, the thrust is set to zero.
    # (The thrust is only important to determine the largest difference between thrust and drag,
    # if the required thrust is higher than the available thrust, the climb cannot be sustained)
    for y in np.arange(0,60,y_step):
        flight_path_angle = y
        flight_path_angle *= Q_("deg")
        L_req = W * np.cos(flight_path_angle)
        alpha_req = L_req / (0.5 * rho * V ** 2 * S * C_L_alpha)
        C_D = C_d_0 + (C_L_alpha * alpha_req) ** 2 / (m.pi * A * e)
        D = C_D * 0.5 * rho * V ** 2 * S
        T_req = D + W * np.sin(flight_path_angle)
        #T = min(min(P_to * eta_prop / V, Tmax).magnitude,T_req.magnitude)
        P_req = T_req * V / eta_prop
        P_req.ito(Q_("hp"))
        P_ratio = P_req / P_to
        alpha_req.ito(Q_("deg"))
        #print(alpha_req)

        if T_req >= (min(P_to * eta_prop / V, Tmax)) or alpha_req > alpha_max:
            T = Q_("0 N")
        else:
            T = T_req

        D_list.append(D.magnitude)
        T_list.append(T.magnitude)
        T_D_list.append((T-D).magnitude)
        y_list.append(y)

    # Gets the largest difference between thrust and drag for the speed concerned
    T_D_max = T_D_list.index(max(T_D_list))
    # print (T_D_max)
    # print("The maximum sustained climb angle is ", y_list [T_D_max])
    y = y_list [T_D_max]
    y_max_list.append(y)
    V_list.append(V.magnitude)

#plt.figure(1)
#plt.plot(V_list, D_list, color = 'r')
#plt.plot(V_list, T_list)
#plt.ylabel('Drag and Trust [N]')
#plt.xlabel('Velocity [m/s] ')

#plt.figure(2)
#plt.title(("Drag and trust for a velocity of", V))
#plt.plot(y_list, D_list, color = 'r')
#plt.plot(y_list, T_list)
#plt.ylabel('Drag and Trust [N]')
#plt.xlabel('Climb angle [deg] ')
#plt.axhline(min(P_to * eta_prop / V, Tmax).magnitude)

plt.figure(3)
plt.plot(V_list, y_max_list, color = 'r')
plt.ylabel('Max sustained climb angle [deg]')
plt.xlabel('Velocity [m/s] ')


plt.show()