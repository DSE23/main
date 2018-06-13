import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
import time
from Misc import ureg, Q_
from Geometry import Geometry as GM
from Aerodynamics import General as Ageneral
from Aerodynamics import Wing as AWing
import Performance as PF
from Misc import Init_parm as IP
import matplotlib.pyplot as plt
import math as m



# Get parameters
P_to = PF.P_to
C_d_0 = Ageneral.CD_0
mass = GM.Masses.W_MTOW
W= mass * Q_("9.81 m/s**2")
S = GM.Wing.S
rho = Q_("1.225 kg/(m**3)")
C_l_max = AWing.CL_max
A = GM.Wing.A
e = AWing.Oswald_e
eta_prop = PF.eta_prop
dp = PF.dp
V_stall =  PF.V_stall_clean
alpha= AWing.alpha_stall
V_lof= 1.05 * V_stall
Tmax = (P_to**2*eta_prop**2*m.pi*dp**2/2*rho)**(1/3)


t_list=[]
y_list=[]
x_list=[]
V_list =[]

p= Q_("2.75 N/m**2")

dt= Q_("0.1 s")
V= Q_("0.0001 m/s")
Vx= Q_("0 m/s")
Vy= Q_("0 m/s")
x=0
t=0


while V < V_lof:
    T = min(P_to * eta_prop / V, Tmax)
    mu = 0.04 #(0.005 + (1 / p) * (0.01 + 0.0095 * (V / 100) ** 2)).magnitude
    C_d = C_d_0
    D = C_d * 0.5 * rho * V ** 2 * S + mu * (W)

    Fx = T-D
    Fy =0

    acceleration_x = Fx/mass
    acceleration_y = 0

    Vx += acceleration_x* dt
    V = Vx
    x += Vx* dt + 0.5 * acceleration_x * dt**2

    t += dt

    V_list.append(Vx.magnitude)
    t_list.append(t.magnitude)
    x_list.append(x)

print("Take-off ground run",x, "time to do this:", t)
plt.plot(t_list, V_list)
plt.ylabel('Velocity  [m/s]')
plt.xlabel('time [s] ')

plt.show()








