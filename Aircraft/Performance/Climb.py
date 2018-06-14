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
alpha = AWing.alpha_stall
V_lof = 1.05 * V_stall
V_max = PF.V_a_clean
Tmax = (P_to**2*eta_prop**2*m.pi*dp**2/2*rho)**(1/3)

flight_path_angle = Q_("45 deg")
flight_path_angle.ito(Q_("rad"))

V_list =[]
alpha_list = []
P_list = []



dV = Q_("1 m/s")

for V in np.arange(V_stall.magnitude, V_max.magnitude, 1):
    V *= Q_("m/s")

    L_req = W * np.cos(flight_path_angle)
    alpha_req = L_req / (0.5 * rho * V**2 * S * C_L_alpha)
    C_D = C_d_0 + (C_L_alpha * alpha_req)**2/ (m.pi * A * e)
    D = C_D * 0.5 * rho * V**2 *S
    T_req = D + W * np.sin(flight_path_angle)
    P_req = T_req * V / eta_prop
    P_ratio = P_req / P_to
    alpha_req.ito(Q_("deg"))

    print(V, alpha_req, P_ratio)