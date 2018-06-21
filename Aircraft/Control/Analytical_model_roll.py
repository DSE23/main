import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path

import numpy as np
from Misc import ureg, Q_      # Imports the unit registry fron the Misc folder
from Geometry import Geometry as GM
from Misc import Init_parm as IP
from Performance import Performance as PF

import math as m
import matplotlib.pyplot as plt


rho = Q_("1.225 kg/m**3")                    # Density at sea level
V_stall = PF.V_stall_clean
V_a = PF.V_a_clean
d_delta_a = GM.Wing.delta_a *2         #rad delta aileron form minus to plus
d_delta_a.ito(Q_("rad"))
d_s_a = Q_("0.41 m")                 #stick deflection
Sa = GM.Wing.S_a                       #m2     aileron surface area
ca = GM.Wing.c_a                       #m     aileron cord
C_h_alpha = IP.C_h_alpha         #hinge moment
C_l_delta_a = IP.C_l_delta_a
C_l_p = IP.C_l_p
C_h_delta = IP.C_h_delta
b = GM.Wing.b
maxdefl_aileron = GM.Wing.delta_a
maxdefl_aileron.ito(Q_("rad"))
Fa= Q_("1000 N")                       #Newton stick force


dt = 0.01
V_list = []
p_list = []
aileron_deflection_list=[]
for V in np.arange(V_stall.magnitude, V_a.magnitude , 0.5):
    V *= Q_("m/s")
    t = 0
    p = Q_("0 rad/s")

    while t < 1:
        delta_alpha_a = (p * 0.5 * b / V) #min(p * 0.5 * b / V, 10)  # delta angle of attack at the tip
        delta_aileron_force = (Fa / (
                    -d_delta_a / d_s_a * 0.5 * rho * V ** 2 * Sa * ca) - C_h_alpha * delta_alpha_a) * 2 / C_h_delta  # aileron deflection for force, hingemoment and speed
        aileron_deflection = min(abs(maxdefl_aileron), abs(delta_aileron_force))
        #print(aileron_deflection)

        p = - C_l_delta_a / C_l_p * aileron_deflection * 2 * V / b
        t = t + dt
    #print (aileron_deflection/m.pi*180)
    #print (delta_alpha_a/m.pi*180)
    p_list.append(p.magnitude)
    V_list.append(V.magnitude)
    aileron_deflection_list.append(aileron_deflection)

plt.figure(1)
plt.plot(V_list, p_list)
plt.ylabel('Roll rate [rad/sec]')
plt.xlabel('Speed [m/s] ')
plt.axhline(y=7.86, color='r')

plt.figure(2)
plt.plot(V_list, aileron_deflection_list)
plt.ylabel('Aileron deflection [rad]')
plt.xlabel('Speed [m/s] ')

plt.show()

