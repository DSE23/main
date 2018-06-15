# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7
Author: Jordy van Leeuwen

This code discretizes the main wing, horizontal tail and vertical tail.
Using non-linear EOM and lookup-tables for the aerodynamic properties, the
manouevre rates and resulting moments are determined.
"""

import sys, os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')  # This suppresses the print commands from the imported modules
sys.path.append('../')              # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
from Geometry import Geometry
from Aerodynamics import Wing as Awing
from Inertia import Inertia
from Misc import Init_parm as IP
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import math
import pandas as pd
import time
sys.stdout = old_stdout             # This enables printing again
t0 = time.time()

# Variables
l_a = Q_("0.3015 m")        # Set aileron length
n_of_disc_w = 30            # number of parts wing is discretized
n_of_disc_h = 10            # number of parts HT is discretized
n_of_disc_v = 10            # number of parts VT is discretized
da = Q_("30 deg")           # aileron deflection
dr = Q_("0 deg")           # rudder deflection
de = Q_("0 deg")           # elevator deflection
alpha_nose = Q_("0.03 rad") # angle of attack of nose
beta_nose  = Q_("0. rad")   # angle of sideslip of nose
V_inf = Q_("110.3 m/s")     # V infinity
t_current = Q_("0 s")       # Start time of sim
dt = Q_("0.01 s")           # Time step of sim
t_sim = Q_("0.5 s")         # End time of sim
l_h = Q_("3.6444 m")        # Tail arm ac-ac
p = Q_("0. 1/s")            # initial roll rate  [rad/s]
q = Q_("0. 1/s")            # initial pitch rate [rad/s]
r = Q_("0. 1/s")            # initial yaw rate   [rad/s]
Phi   = Q_("0. rad")        # Initial euler angle around x-axis
Psi   = Q_("0. rad")        # Initial euler angle around z-axis
Theta = Q_("0. rad")        # Initial euler angle around y-axis

# Definitions
def local_chord(z, c_r, c_t, half_b):
    # Calculates the chord at location z(distance from center)
    return c_r - (c_r - c_t) / half_b * z


# import airfoil lookup tables
data = pd.read_csv('aerodynamic_data_ms15.dat', ' ', header=None).values

def lookup_data(alpha, ca_c, da, chord):
    # Looksup and interpolates Cl and Cd based on alpha, ca_c and da
    alpha = math.degrees(alpha)
    indexda = abs(da)//5
    indexca_c = int(100*ca_c)-1
    localdata = data[int(indexda*50*51+indexca_c*51):int(indexda*50*51+(indexca_c+1)*51),:]
    non_zero_max = max(np.argwhere(localdata[:, 0]))[0]  # last non-zero row
    localdata = localdata[:non_zero_max+1,:]
    Cl_local = interpolate.interp1d(localdata[:,0], localdata[:,1], 'linear', fill_value='extrapolate')
    Cd_local = interpolate.interp1d(localdata[:,0], localdata[:,2], 'linear', fill_value='extrapolate')
    Cm_local = interpolate.interp1d(localdata[:,0], localdata[:,4], 'linear', fill_value='extrapolate')
    
    if da > 0:
        Cl = Cl_local(alpha)
        Cd = Cd_local(alpha)
        Cm = Cm_local(alpha)
    else:
        Cl = -Cl_local(-alpha)
        Cd = Cd_local(-alpha)
        Cm = Cm_local(-alpha)
    
    Cn = Cl*math.cos(alpha) + Cd*math.sin(alpha)
    if Cn == 0.:
        xcp = 0.
    else:
        xcp = 0.25-Cm/Cn
    return Cl, Cd, xcp

# Import other forces
T = 0.*ureg.N
D_fus_gear = ((0.01998*V_inf*V_inf).magnitude)*ureg.N
W = (IP.MTOW*9.80665)*ureg.N

# Import Aircraft Geometry
I_yy = Inertia.I_yy
I_xx = Inertia.I_xx
I_zz = Inertia.I_zz
b_w = Geometry.Wing.b
S_w = Geometry.Wing.S
c_r_w = Geometry.Wing.c_r
c_t_w = Geometry.Wing.c_t
AR_w = Geometry.Wing.A
e_w = Awing.Oswald_e
b_h = Geometry.H_tail.b
S_h = Geometry.H_tail.S
c_r_h = Geometry.H_tail.c_r
c_t_h = Geometry.H_tail.c_t
AR_h = Geometry.H_tail.A
b_v = Geometry.V_tail.b
S_v = Geometry.V_tail.S
c_r_v = Geometry.V_tail.c_r
c_t_v = Geometry.V_tail.c_t
AR_v = Geometry.V_tail.A
t_c_v = Geometry.V_tail.t_c
Z_v = Geometry.V_tail.Z_v

cabin_width = Geometry.Fuselage.cabin_w + 0.1 * ureg.m  # import cabin width
vt_width = t_c_v * c_r_v + 0.1 * ureg.m                 # import width of VT
rho = IP.rho0                                           # import density

bloc_w = (b_w-cabin_width) / n_of_disc_w  # span of each station wing
bloc_h = b_h / n_of_disc_h  # Span of each station HT
bloc_v = b_v / n_of_disc_v  # Span of each station VT
half_b_w = b_w / 2          # half span wing
half_b_h = b_h / 2          # half span HT

V_s = IP.V_stall_clean              # stall speed
V_a = IP.V_a_clean                  # manoeuvring speed

# Setup lists with station boundaries
n_chords_w = int(n_of_disc_w / 2)
n_chords_h = int(n_of_disc_h / 2)
kwlst = n_chords_w * [0] + [-cabin_width / 2, 0 * ureg.m, cabin_width / 2] + n_chords_w * [0]
khlst = n_chords_h * [0] + [-vt_width / 2, 0 * ureg.m, vt_width / 2] + n_chords_h * [0]
kvlst = (n_of_disc_v + 1) * [0]
for j in range(n_chords_w):
    kwlst[j] = bloc_w * j - half_b_w
    kwlst[-1 - j] = half_b_w - bloc_w * j
for l in range(n_chords_h):
    khlst[l] = bloc_h * l - half_b_h
    khlst[-1 - l] = half_b_h - bloc_h * l
for k in range(n_of_disc_v + 1):
    kvlst[k] = Z_v - b_v + bloc_v * k
    kvlst[k] = Z_v - b_v + bloc_v * k


# Init sim
# n_V = 0
#Vrange = np.arange(V_s, V_a, 1)
#pmax = np.transpose(np.vstack([Vrange, np.zeros((1, len(Vrange)))[0]]))
# for V in Vrange:
#    V_inf = V*ureg.m/ureg.s

running = True
n = 0
plst  = np.zeros((1, int((t_sim / dt).magnitude)))[0]
pdlst = np.zeros((1, int((t_sim / dt).magnitude)))[0]
Rlst  = np.zeros((1, int((t_sim / dt).magnitude)))[0]
tlst  = np.arange(0, t_sim.magnitude, dt.magnitude)

while running:
    if (t_current.magnitude + dt.magnitude) > t_sim.magnitude:
        running = False
    
    disc_wing_w = np.zeros((len(kwlst)-1, 8))  # 2D array discretized wing
    disc_wing_h = np.zeros((len(khlst)-1, 6))  # 2D array discretized HT
    disc_wing_v = np.zeros((len(kvlst)-1, 6))  # 2D array discretized VT

    disc_wing_w[(range(n_chords_w)), 0] = da  # set aileron for postive stations
    disc_wing_w[(range(int(n_of_disc_w + 2 - n_chords_w), n_of_disc_w + 2)), 0] = -da  # set aileron for negative stations
    disc_wing_h[(range(n_chords_h)),0] = de
    disc_wing_h[(range(int(n_of_disc_h + 2 - n_chords_h), n_of_disc_h + 2)), 0] = -de
    disc_wing_v[:,0] = dr
    for i in range(0, len(kwlst)-1):
        # calculate lift and drag for discretized wing
        da_local = disc_wing_w[i][0]
        b1 = kwlst[i]                                       # Y boundary left of piece
        b2 = kwlst[i+1]                                     # Y boundary right of piece
        y_i = (b1 + b2) / 2                                 # Y centre of piece
        c1 = local_chord(abs(b1), c_r_w, c_t_w, half_b_w)   # Chord left of piece
        c2 = local_chord(abs(b2), c_r_w, c_t_w, half_b_w)   # Chord right of piece
        cc = local_chord(abs(y_i), c_r_w, c_t_w, half_b_w)  # Chord centre of piece
        ca_c = (l_a / cc).magnitude                         # percentage aileron chord over local aileron
        Sloc = (c1 + c2) / 2 * (b2 - b1)                    # Surface area of piece

        alpha_w = alpha_nose + p * y_i / V_inf
        alpha_i = 0.
        running_alpha_i = True
        while running_alpha_i:
            alpha_e = alpha_w - alpha_i
            Cl, Cd, Cp = lookup_data(alpha_e, ca_c, da_local,cc)
            alpha_i_new = Cl / (math.pi * AR_w * e_w)
            if (alpha_i_new-alpha_i)/alpha_i_new < 0.01:
                running_alpha_i = False
            alpha_i = alpha_i_new
        alpha = alpha_w - alpha_i
        downwash_angle = 2 * Cl / (math.pi * AR_w)
        Cl, Cd, Cp = lookup_data(alpha, ca_c, da_local,cc)
        Cdi = Cl * Cl / (math.pi * AR_w * e_w)
        Cd = Cd + Cdi
        Lift = 0.5 * rho * V_inf ** 2 * Sloc * Cl
        Drag = 0.5 * rho * V_inf ** 2 * Sloc * Cd
        disc_wing_w[i][6] = Cp.magnitude
        disc_wing_w[i][1] = Sloc.magnitude
        disc_wing_w[i][2] = Drag.magnitude
        disc_wing_w[i][3] = Lift.magnitude
        disc_wing_w[i][4] = (Drag * y_i).magnitude
        disc_wing_w[i][5] = (Lift * y_i).magnitude

    for i in range(0, len(khlst)-1):
        # calculate lift and drag for discretized HT
        de_local = disc_wing_h[i][0]
        b1 = khlst[i]
        b2 = khlst[i+1]
        y_i = (b1 + b2) / 2
        c1 = local_chord(abs(b1), c_r_h, c_t_h, half_b_h)
        c2 = local_chord(abs(b2), c_r_h, c_t_h, half_b_h)
        cc = local_chord(abs(y_i),c_r_h, c_t_h, half_b_h)
        Sloc = (c1 + c2) / 2 * (b2 - b1)
        alpha = alpha_nose + p*y_i/V_inf
        Cl, Cd, Cp = lookup_data(alpha, 0.01, de_local,cc)
        Cdi = Cl * Cl / (math.pi * AR_w * e_w)
        #Cd = Cd + Cdi
        Lift = 0.5 * rho * V_inf ** 2 * Sloc * Cl
        Drag = 0.5 * rho * V_inf ** 2 * Sloc * Cd
        disc_wing_h[i][1] = Sloc.magnitude
        disc_wing_h[i][2] = Drag.magnitude
        disc_wing_h[i][3] = Lift.magnitude
        disc_wing_h[i][4] = (Drag * y_i).magnitude
        disc_wing_h[i][5] = (Lift.magnitude * y_i.magnitude)

    for i in range(0, len(kvlst)-1):
        # calculate lift and drag for discretized VT
        dr_local = disc_wing_v[i][0]
        b1 = kvlst[i]
        b2 = kvlst[i+1]
        z_i = (b1 + b2) / 2
        c1 = local_chord(abs(b1 - Z_v), c_r_v, c_t_v, b_v)
        c2 = local_chord(abs(b2 - Z_v), c_r_v, c_t_v, b_v)
        Sloc = (c1 + c2) / 2 * (b2 - b1)
        beta = beta_nose + p * z_i / V_inf
        Cl, Cd, Cp = lookup_data(beta, 0.01, dr_local,1.0*ureg.m)
        Lift = 0.5 * rho * V_inf ** 2 * Sloc * Cl
        Drag = 0.5 * rho * V_inf ** 2 * Sloc * Cd
        disc_wing_v[i][1] = Sloc.magnitude
        disc_wing_v[i][2] = beta.magnitude
        disc_wing_v[i][3] = Lift.magnitude
        disc_wing_v[i][4] = c2.magnitude
        disc_wing_v[i][5] = (Lift.magnitude * z_i.magnitude)
        
    Sum_Dw_y = sum(disc_wing_w[:,4])*ureg.N*ureg.m
    Sum_Dh_y = sum(disc_wing_h[:,4])*ureg.N*ureg.m
    Sum_Lw_y = sum(disc_wing_w[:,5])*ureg.N*ureg.m
    Sum_Lh_y = sum(disc_wing_h[:,5])*ureg.N*ureg.m
    Sum_Lv_y = sum(disc_wing_v[:,5])*ureg.N*ureg.m
    Sum_Lw   = sum(disc_wing_w[:,3])*ureg.N
    Sum_Lh   = sum(disc_wing_h[:,3])*ureg.N
    Sum_Dw   = sum(disc_wing_w[:,2])*ureg.N
    Sum_Dh   = sum(disc_wing_h[:,2])*ureg.N
    Sideforce= sum(disc_wing_v[:,3])*ureg.N
    
    L = - Sum_Lw_y - Sum_Lh_y - Sum_Lv_y 
    R = (Sum_Dw_y + Sum_Dh_y)*math.cos(beta_nose) + Sideforce*l_h   # INCOMPLETE!
    M = - Sum_Lh * l_h*math.cos(alpha_nose)                         # INCOMPLETE!
    Fx = T - math.cos(beta_nose)*(Sum_Dw + Sum_Dh + D_fus_gear) - W*math.sin(Theta)
    Fy = Sideforce - math.sin(beta_nose)*(Sum_Dw + Sum_Dh + D_fus_gear) + W*math.sin(Theta)
    Fz = - Sum_Lw - Sum_Lh + W*math.cos(Phi)
    
    p_dot = L / I_xx
    p    += p_dot * dt
    Phi  += p *dt
    
    q_dot  = M / I_yy
    q     += q_dot * dt
    Theta += q * dt
    
    r_dot = R / I_zz
    r    += r_dot * dt
    Psi  += r * dt
    
    # update lists for plots
    plst[n]  = p.magnitude
    pdlst[n] = p_dot.magnitude
    Rlst[n]  = R.magnitude

    t_current += dt
    n += 1

#    pmax[n_V,1] = max(plst)
#    n_V += 1

# plt.plot(pmax[:,0],pmax[:,1])
# plt.plot(tlst,plst)
plt.plot(tlst, np.degrees(plst))
plt.plot(tlst, pdlst)
#plt.show()

print("Finished in:", round(time.time() - t0, 1), "s")