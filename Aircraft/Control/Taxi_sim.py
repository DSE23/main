# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:03:58 2018

@author: jurian
"""


import sys
sys.path.append('../')   # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
from matplotlib import pyplot as plt
import math as m
import numpy as np
from Geometry import Geometry
from Aerodynamics import Aeroprops
from Performance import Performance

# First check CG position for tipping:

Z_mainlg = Geometry.Landing_gear.Z_mainlg
X_mainlg = Geometry.Landing_gear.X_mainlg
X_taillg = Geometry.Landing_gear.X_taillg
cgangle_fw = Q_("16 deg")
cgangle_rear = Q_("25 deg")
X_maxfront = X_mainlg + Z_mainlg * np.tan(cgangle_fw)
X_maxback = X_mainlg + Z_mainlg * np.tan(cgangle_rear)
X_cgmtow = Geometry.CG.CG_mtow
X_cgoew = Geometry.CG.CG_OEW
V_taxi = Q_("5 kts")
V_taxi.ito(Q_("m/s"))
rho0 = Q_("1.225 kg/m**3")
g0 = Performance.g0.magnitude
g0 = g0 * Q_("m/s**2")
m_mtow = Geometry.Masses.W_MTOW
W_mtow =  m_mtow * g0
V_taxi = Q_("5 kts")
V_taxi.ito(Q_("m/s"))
rho0 = Q_("1.225 kg/m**3")
g0 = Performance.g0.magnitude
g0 = g0 * Q_("m/s**2")
m_mtow = Geometry.Masses.W_MTOW
W_mtow = m_mtow * g0
Wheel_base = X_mainlg - X_taillg
if not X_maxfront < X_cgmtow < X_maxback or not X_maxfront < X_cgoew < X_maxback:
    print('\x1b[3;37;41m' + "Front landing gear does not meet requirements" + '\x1b[0m')

Stat_load_main = (-W_mtow * (X_taillg - X_cgmtow)/Wheel_base)/2
Stat_load_tail = (W_mtow * (X_mainlg - X_cgmtow)/Wheel_base)


# Calculate turn rate at 5 kts

alpha = np.radians(Geometry.Landing_gear.Tip_angle)
CL_alphaw = Aeroprops.CL_alpha_wing
CL_alphah = Aeroprops.CL_alpha_ht
S_wing = Geometry.Wing.S
S_htail = Geometry.H_tail.S
L_wing = (CL_alphaw * alpha) * 0.5 * rho0 * V_taxi**2 * S_wing
L_htail = (CL_alphah * alpha) * 0.5 * rho0 * V_taxi**2 * S_htail
XWing_25c = Geometry.CG.X_wing
XHtail_25c = Geometry.H_tail.X_h
L_tot = L_wing + L_htail
N_force = W_mtow  - L_tot
M_lift = - L_wing * (XWing_25c - X_cgmtow) - L_htail * (XHtail_25c - X_cgmtow)
N_main = -N_force * (X_taillg - X_cgmtow)/Wheel_base - M_lift/Wheel_base
N_tail = N_force * (X_mainlg - X_cgmtow)/Wheel_base + M_lift/Wheel_base
Max_allow_p = Q_("35 psi") / g0
Max_allow_p.ito("kg/cm**2")
print(Max_allow_p)
print(Stat_load_main / g0)
print(Stat_load_tail / g0)
t_end = Q_("10 s")
dt = Q_("0.01 s")
Psy = Q_("0 deg")                             # Initial heading
Y = 0                                       # Initial lateral position
X = 0                                       # Initial X-position
delta_w = Q_("0 deg")                       # Initial Wheel deflection
Psilist =[]
Psidotlist =[]
Xlist = []
Ylist = []
tlist =[]
for t in np.arange(0, t_end.magnitude + dt.magnitude, dt.magnitude):
    if t >= 8:
        delta_w = Q_("0 deg")
    elif t >= 2:
        delta_w = Q_("25 deg")
    Psydot = V_taxi/Wheel_base * np.tan(delta_w)    
    Ydot = V_taxi * np.sin(Psy)
    Xdot = V_taxi * np.cos(Psy)
    Psy = Psy + Psydot*dt
    Y = Y + Ydot*dt
    X = X + Xdot*dt
    Psilist.append(Psy.magnitude)
    Psidotlist.append(Psydot)
    Xlist.append(X.magnitude)
    Ylist.append(Y.magnitude)
    tlist.append(t)
plt.plot(tlist, Psilist)
print(min(Psidotlist))
plt.show()

    
    


