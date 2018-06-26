import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
from Misc import ureg, Q_
from Geometry import Geometry as GM
from Aerodynamics import Aeroprops as Aeroprops
from Aerodynamics import Wing as AWing
from Performance import Performance as PF
from Propulsion_and_systems import Propeller as Prop
from Misc import Init_parm as IP
import matplotlib.pyplot as plt
import math as m



# Get parameters
P_to = PF.P_to.magnitude
P_to = P_to * Q_("kg*m**2/s**3")
C_d_0 = Aeroprops.CD0_tot
mass = GM.Masses.W_MTOW
W = mass * Q_("9.81 m/s**2")
S = GM.Wing.S
rho = Q_("1.225 kg/(m**3)")
C_L_alpha = AWing.CL_alpha
C_L_alpha *= Q_("1/rad")
C_l_max = AWing.CL_max_hld
A = GM.Wing.A
e = AWing.Oswald_e
eta_prop = PF.eta_prop
dp = PF.dp
V_stall =  PF.V_stall_hld.magnitude
V_stall *=Q_("m/s")
alpha_max = AWing.alpha_stall
alpha_max *= Q_("deg")
V_REF = 1.3 * V_stall           # from CS23.73
Tmax = Prop.Tstatic.magnitude 
Tmax *= Q_("N") 

g = PF.g0.magnitude
g *= Q_("m/s**2")

# initial conditions
V = V_REF
h = Q_("15 m")      # from CS23.73, screen height
gradient_of_descent = Q_("3 deg")
x = Q_("0 m")
t = Q_("0 s")


# inputs
n = 1.1     # Load factor at landing
dt = Q_("0.01 s")
mu = 0.4

# Empty lists
t_list =[]
h_list =[]
x_list =[]
V_list =[]

# setting needed for desired glide slope
flight_path_angle = - gradient_of_descent
L_req = W * np.cos(flight_path_angle)
alpha_req = L_req / (0.5 * rho * V**2 * S * C_L_alpha)
C_D = C_d_0 + (C_L_alpha * alpha_req)**2/ (m.pi * A * e)
D = C_D * 0.5 * rho * V**2 *S
T_req = D + W * np.sin(flight_path_angle)
P_req = T_req * V / eta_prop
P_req.ito(Q_("hp"))
P_ratio = P_req / P_to
alpha_req.ito(Q_("deg"))
R = V**2 /(g*(n-1))

h_f = R - R * np.cos(flight_path_angle)
print( "This program runs for a load factor of ", n, " at landing with a flight path angle of ", flight_path_angle)
print("Flare happens at ", round(h_f,2), ". The required AoA is " ,  round(alpha_req), "with a required power of", round(P_req,2))

while h > h_f:
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
    #flight_path_angle_dot = Fy / (mass * V)

    x += V * np.cos(flight_path_angle) * dt
    h += V * np.sin(flight_path_angle) * dt

    V += V_dot * dt
    #flight_path_angle += flight_path_angle_dot * dt
    #print(flight_path_angle)

    t += dt
    V_list.append(V.magnitude)
    t_list.append(t.magnitude)
    h_list.append(h.magnitude)
    x_list.append(x.magnitude)

t2 = (2 * m.pi * R / V) / Q_("360 deg") * flight_path_angle
flight_path_angle_dot = flight_path_angle / t2


print("The change of flight path angle needed for the flare is " , round(flight_path_angle_dot,2), ". The flare takes ", round(-t2,2))


while flight_path_angle < Q_("0 deg"):

    flight_path_angle = flight_path_angle + flight_path_angle_dot * dt
    x += V * np.cos(flight_path_angle) * dt
    h += V * np.sin(flight_path_angle) * dt

    t += dt
    V_list.append(V.magnitude)
    t_list.append(t.magnitude)
    h_list.append(h.magnitude)
    x_list.append(x.magnitude)
print(V)
print("The distance covered from screen height to touchdown is", round(x,2))

while V > Q_("1 m/s"):
    T = 0 # min(P_to * eta_prop / V, Tmax)
    C_d = C_d_0
    D = C_d * 0.5 * rho * V ** 2 * S + mu * (W)

    Fx = T-D
    Fy = 0

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


print("The total landing distance is",x)
print("The total landing time is",t)
print(V_REF)

plt.figure(1)
plt.plot(t_list, h_list)
plt.ylabel('height [m]')
plt.xlabel('time [s] ')

plt.figure(2)
plt.plot(x_list, h_list)
plt.ylabel('Height [m]')
plt.xlabel('Distance [m] ')
plt.axis('scaled')
plt.xlim(0,x.magnitude)
plt.ylim(0,25)

plt.figure(3)
plt.plot(x_list, V_list)
plt.ylabel('velocity [m/s]')
plt.xlabel('distance [m] ')

#plt.xlim(200,300)
#plt.ylim(0,100)

plt.show()