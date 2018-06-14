import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
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
V_lof = 1.05 * V_stall
Tmax = (P_to**2*eta_prop**2*m.pi*dp**2/2*rho)**(1/3)

#Set parameters
t_run = Q_("100 s")
dt = Q_("0.01 s")
alpha_input = Q_("2 deg")

t_list =[]
h_list =[]
x_list =[]
V_list =[]

p = Q_("2.75 N/m**2")

# initial conditions
V = Q_("0.0001 m/s")
x = Q_("0 m")
h = Q_("0 m")
t = Q_("0 s")


while V < V_lof:
    T = min(P_to * eta_prop / V, Tmax)
    mu = 0.04 #(0.005 + (1 / p) * (0.01 + 0.0095 * (V / 100) ** 2)).magnitude
    C_d = C_d_0
    D = C_d * 0.5 * rho * V ** 2 * S + mu * (W)

    Fx = T-D
    Fy =0

    acceleration_x = Fx/mass
    acceleration_y = 0

    V += acceleration_x* dt
    x += V* dt + 0.5 * acceleration_x * dt**2
    h = Q_("0 m ")

    t += dt
    V_list.append(V.magnitude)
    t_list.append(t.magnitude)
    h_list.append(h.magnitude)
    x_list.append(x.magnitude)


P = P_to
flight_path_angle = Q_("0 deg")
flight_path_angle.ito(Q_("rad"))

print(h)
print("Take-off ground run",x, "time to do this:", t)


while flight_path_angle < Q_(" 45 deg "):
    alpha = alpha_max
    alpha.ito(Q_("rad"))
    C_L = C_L_alpha * alpha
    L = 0.5 * rho * V**2 *S * C_L
    C_d = C_d_0 + C_L**2 / (m.pi*A*e)
    D = C_d * 0.5 * rho * V**2 * S
    T = min(P * eta_prop / V, Tmax)

    Fx = T - D - W * np.sin(flight_path_angle)
    Fy = L - W * np.cos(flight_path_angle)


    V_dot = Fx / mass
    flight_path_angle_dot = Fy / (mass * V)

    x += V * np.cos(flight_path_angle) * dt
    h += V * np.sin(flight_path_angle) * dt

    V += V_dot * dt
    flight_path_angle += flight_path_angle_dot * dt
    #print(flight_path_angle)

    t += dt
    V_list.append(V.magnitude)
    t_list.append(t.magnitude)
    h_list.append(h.magnitude)
    x_list.append(x.magnitude)

L_req = W * np.cos(flight_path_angle)
alpha_req = L_req / (0.5 * rho * V**2 * S * C_L_alpha)
C_D = C_d_0 + (C_L_alpha * alpha_req)**2/ (m.pi * A * e)
D = C_D * 0.5 * rho * V**2 *S
T_req = D + W * np.sin(flight_path_angle)
P_req = T_req * V / eta_prop
P_ratio = P_req / P_to
alpha_req.ito(Q_("deg"))

print(alpha_req, P_ratio)


while t < Q_(" 60 s "):
    alpha = alpha_req
    P = P_req
    alpha.ito(Q_("rad"))
    C_L = C_L_alpha * alpha
    L = 0.5 * rho * V**2 *S * C_L
    C_d = C_d_0 + C_L**2 / (m.pi*A*e)
    D = C_d * 0.5 * rho * V**2 * S
    T = min(P * eta_prop / V, Tmax)

    Fx = T - D - W * np.sin(flight_path_angle)
    Fy = L - W * np.cos(flight_path_angle)


    V_dot = Fx / mass
    flight_path_angle_dot = Fy / (mass * V)

    x += V * np.cos(flight_path_angle) * dt
    h += V * np.sin(flight_path_angle) * dt

    V += V_dot * dt
    flight_path_angle += flight_path_angle_dot * dt
    #print(flight_path_angle)

    t += dt
    V_list.append(V.magnitude)
    t_list.append(t.magnitude)
    h_list.append(h.magnitude)
    x_list.append(x.magnitude)


plt.figure(1)
plt.plot(t_list, V_list)
plt.ylabel('Velocity  [m/s]')
plt.xlabel('time [s] ')

plt.figure(2)
plt.plot(t_list, h_list)
plt.ylabel(' height [m]')
plt.xlabel('time [s] ')

plt.figure(3)
plt.plot(x_list, h_list)
plt.ylabel(' height [m]')
plt.xlim(0,250)
plt.ylim(0,250)
plt.xlabel('distance [m] ')

plt.show()