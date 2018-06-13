import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path

import numpy as np
from Misc import ureg, Q_      # Imports the unit registry fron the Misc folder
from Geometry import Wing
from Misc import Init_parm as IP

import math as m
import matplotlib.pyplot as plt


rho0 = 1.225                    # Density at sea level
V_stall = IP.V_stall
V_a = IP.V_man
d_delta_a= IP.d_delta_a         #rad delta aileron form minus to plus
d_s_a= IP.d_s_a                 #stick deflection
Sa= IP.Sa                       #m2     aileron surface area
ca= IP.ca                       #m     aileron cord
C_h_alpha= IP.C_h_alpha         #hinge moment
C_l_delta_a= IP.C_l_delta_a
C_l_p= IP.C_l_p
C_h_delta= IP.C_h_delta
b= IP.b
maxdefl_aileron = IP.maxdefl_aileron #rad
alpha_max= IP.alpha_max         #deg
Fa= IP.Fa                       #Newton stick force
rho= IP.rho0

dt = 0.01
V_list = []
p_list = []

for V in np.arange(V_stall, V_a, 0.5):
    t = 0
    p = 0

    while t < 1:
        delta_alpha_a = (p * 0.5 * b / V) #min(p * 0.5 * b / V, 10)  # delta angle of attack at the tip
        delta_aileron_force = (Fa / (
                    -d_delta_a / d_s_a * 0.5 * rho * V ** 2 * Sa * ca) - C_h_alpha * delta_alpha_a) * 2 / C_h_delta  # aileron deflection for force, hingemoment and speed
        aileron_deflection = min(abs(maxdefl_aileron), abs(delta_aileron_force))

        p = - C_l_delta_a / C_l_p * aileron_deflection * 2 * V / b
        t = t + dt
    #print (aileron_deflection/m.pi*180)
    #print (delta_alpha_a/m.pi*180)
    p_list.append(p)
    V_list.append(V)

plt.figure(1)
plt.plot(V_list, p_list)
plt.ylabel('Roll rate [rad/sec]')
plt.xlabel('Speed [m/s] ')
plt.axhline(y=7.86, color='r')

plt.show()

