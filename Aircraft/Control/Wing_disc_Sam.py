# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7
Author: Jordy van Leeuwen

This code discretizes the main wing, horizontal tail and vertical tail.
Using non-linear EOM and lookup-tables for the aerodynamic properties, the
manouevre rates and resulting moments are determined.

To be updated:
- Trimming the initial condition
- Thrust
- Drag from fuselage and gear
"""

import sys, os
sys.path.append('../')              # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
from Geometry import Geometry
from Aerodynamics import Wing as Awing
from Inertia import Inertia
from Performance import Performance
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import math as m
import time
import pandas as pd
t0 = time.time()

# Variables
l_a = Q_("0.3015 m")        # Set aileron length
l_e = Q_("0.2    m")        # Set elevator length
cr_c= Q_("0.2     ")
n_of_disc_w = 30            # number of parts wing is discretized
n_of_disc_h = 10            # number of parts HT is discretized
n_of_disc_v = 10            # number of parts VT is discretized
da = Q_("0 deg")            # aileron deflection
dr = Q_("0 deg")            # rudder deflection
de = Q_("0 deg")            # elevator deflection
alpha_nose = Q_("0.0195 rad") # angle of attack of nose
beta_nose  = Q_("0. rad")   # angle of sideslip of nose
V_inf = Q_("110.3 m/s")     # V infinity
t_current = Q_("0.0 s")       # Start time of sim
dt = Q_("0.01 s")           # Time step of sim
t_end = Q_("0.1 s")         # End time of sim
l_h = Q_("3.6444 m")        # Tail arm ac-ac horizontal
l_v = Q_("3.7 m")           # Tail arm ac-ac vertical
p = Q_("0. 1/s")            # initial roll rate  [rad/s]
q = Q_("0. 1/s")            # initial pitch rate [rad/s]
r = Q_("0. 1/s")            # initial yaw rate   [rad/s]
Phi   = Q_("0. rad")        # Initial euler angle around x-axis
Psi   = Q_("0. rad")        # Initial euler angle around z-axis
Theta = Q_("0. rad")        # Initial euler angle around y-axis
w = Q_("0. m/s")
u = V_inf
v = Q_("0. m/s")
gamma = Q_("0. rad")
gamma_lateral = Q_("0 rad")
mtow = Geometry.Masses.W_MTOW
g0 = Performance.g0.magnitude * Q_("m/s**2")

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
e_h = e_w
e_v = e_w
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

xcg    = Geometry.CG.CG_mtow
xlemac = Geometry.CG.XLEMAC
MAC = Geometry.Wing.MAC
sweep_LE = Geometry.Wing.Sweep_LE.magnitude*1

cabin_width = Geometry.Fuselage.cabin_w + 0.1 * ureg.m  # import cabin width
vt_width = t_c_v * c_r_v + 0.1 * ureg.m                 # import width of VT
rho = Performance.rho_0.magnitude * Q_("kg/m**3")                        # import density

bloc_w = (b_w-cabin_width) / n_of_disc_w  # span of each station wing
bloc_h = b_h / n_of_disc_h  # Span of each station HT
bloc_v = b_v / n_of_disc_v  # Span of each station VT
half_b_w = b_w / 2          # half span wing
half_b_h = b_h / 2          # half span HT
y_mac = half_b_w*(c_r_w-MAC)/(c_r_w-c_t_w)

V_s = Performance.V_stall_clean            # stall speed
V_a = Performance.V_a_clean                  # manoeuvring speed

# Definitions
def local_chord(z, c_r, c_t, half_b):
    # Calculates the chord at location z(distance from center)
    return c_r - (c_r - c_t) / half_b * z
def trimming(u,ca_c,da, chord):
    # Jurians dclde
    Cl1, dummy, dummy, dummy =  (lookup_data(0.,0.5,0,1.))
    Cl2, dummy, dummy, dummy =  (lookup_data(0.,0.5,10,1.))
    dcl_de = (Cl2-Cl1)/(10)
    
    
    trimming_alpha = True
    alpha_min = m.radians(-15)
    alpha_max = m.radians(15)
    while trimming_alpha:
        alpha_t = (alpha_min+alpha_max)/2.
        

        w = u * m.tan(alpha_t)
        Cl_w,Cd_w,Cm_w,xcp_w = lookup_data(alpha_t, ca_c, da, chord)
        Cl_h,Cd_h,Cm_h,xcp_h = lookup_data(alpha_t, ca_c, de, chord)
        Lift_w = 0.5 * rho * (u*u+w*w) * Cl_w * S_w
        Lift_h = 0.5 * rho * (u*u+w*w) * Cl_h * S_h
        Drag_w = 0.5 * rho * (u*u+w*w) * Cd_w * S_w
        Drag_h = 0.5 * rho * (u*u+w*w) * Cd_h * S_h
        Lift = Lift_w + Lift_h
        Drag = Drag_w + Drag_h
        zforce = -Lift*m.cos(alpha_t) - Drag*m.sin(alpha_t) + W*m.cos(Theta)

        if zforce.magnitude >0.:
            alpha_min = alpha_t
        else:
            alpha_max = alpha_t
        if abs(zforce.magnitude)<1.0:
            trimming_alpha = False
    print (Lift)
    return alpha_t

# import airfoil lookup tables
data = pd.read_csv('aerodynamic_data_ms15.dat', ' ', header=None).values

def lookup_data(alpha, ca_c, da, chord):
    # Looksup and interpolates Cl and Cd based on alpha, ca_c and da
    alpha = m.degrees(alpha)
    indexca_c = int(100*ca_c)-1
    if da % 5 ==0:
        indexda = abs(da)//5
        localdata = data[int(indexda*50*51+indexca_c*51):int(indexda*50*51+(indexca_c+1)*51),:]
    else:
        index1da = abs(da)//5
        index2da = abs(da)//5 + 1
        localdata1 = data[int(index1da*50*51+indexca_c*51):int(index1da*50*51+(indexca_c+1)*51),:]
        localdata2 = data[int(index2da*50*51+indexca_c*51):int(index2da*50*51+(indexca_c+1)*51),:]
        localdata = (localdata2 - localdata1)/5 * da
    non_zero_max = max(np.argwhere(localdata[:, 0]))[0]  # last non-zero row
    localdata = localdata[:non_zero_max+1,:]
    Cl_local = interpolate.interp1d(localdata[:,0], localdata[:,1], 'linear', fill_value='extrapolate')
    Cd_local = interpolate.interp1d(localdata[:,0], localdata[:,2], 'linear', fill_value='extrapolate')
    Cm_local = interpolate.interp1d(localdata[:,0], localdata[:,4], 'linear', fill_value='extrapolate')
    #plt.plot(localdata[:,0],localdata[:,4])
    #plt.plot(localdata[:,0],localdata[:,1])
    
    #plt.plot(localdata[range(1,len(localdata[:,0])),0],0.25-localdata[range(1,len(localdata[:,0])),4]/localdata[range(1,len(localdata[:,0])),1])
    if da >= 0:
        Cl = Cl_local(alpha)
        Cd = Cd_local(alpha)
        Cm = Cm_local(alpha)
    else:
        Cl = -Cl_local(-alpha)
        Cd = Cd_local(-alpha)
        Cm = -Cm_local(-alpha)
    
    Cn = m.cos(m.radians(alpha))*Cl + m.sin(m.radians(alpha))*Cd
    if Cn !=0:
        xcp = 0.25-Cm/Cn
    else:
        xcp = 0.25
    return Cl, Cd, Cm, xcp

# Import other forces
T = Q_("1000 N")
D_fus_gear = ((0.01998*V_inf*V_inf).magnitude)*Q_("N")
W = mtow*g0


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
# empty arrays for plotting in the end:
plst  = np.zeros((1, int((t_end / dt).magnitude)))[0]
pdlst = np.zeros((1, int((t_end / dt).magnitude)))[0]
qlst  = np.zeros((1, int((t_end / dt).magnitude)))[0]
qdlst = np.zeros((1, int((t_end / dt).magnitude)))[0]
Fzlst = np.zeros((1, int((t_end / dt).magnitude)))[0]
tlst  = np.arange(0, t_end.magnitude, dt.magnitude)

#alpha_nose = trimming(V_inf,0.1,0,1)
print ("Alpha_t:",alpha_nose)

for t_current in np.arange(0,(t_end).magnitude,dt.magnitude):
    print(n)
    disc_wing_w = np.zeros((len(kwlst)-1, 20))  # 2D array discretized wing
    disc_wing_h = np.zeros((len(khlst)-1, 20))  # 2D array discretized HT
    disc_wing_v = np.zeros((len(kvlst)-1, 20))  # 2D array discretized VT

    # Setting the control surface deflection for each station
    disc_wing_w[(range(n_chords_w)), 0] = da  
    disc_wing_w[(range(int(n_of_disc_w + 2 - n_chords_w), n_of_disc_w + 2)), 0] = -da  # set aileron for negative stations
    disc_wing_h[(range(n_chords_h)),0] = de
    disc_wing_h[(range(int(n_of_disc_h + 2 - n_chords_h), n_of_disc_h + 2)), 0] = de
    disc_wing_v[:,0] = dr
    
    # Calculate for WING
    for i in range(0, len(kwlst)-1):
        da_local = disc_wing_w[i][0]                        # Local aileron deflection of piece
        b1 = kwlst[i]                                       # Y boundary left of piece
        b2 = kwlst[i+1]                                     # Y boundary right of piece
        y_i = (b1 + b2) / 2                                 # Y centre of piece
        c1 = local_chord(abs(b1), c_r_w, c_t_w, half_b_w)   # Chord left of piece
        c2 = local_chord(abs(b2), c_r_w, c_t_w, half_b_w)   # Chord right of piece
        cc = local_chord(abs(y_i), c_r_w, c_t_w, half_b_w)  # Chord centre of piece
        ca_c = (l_a / cc).magnitude                         # percentage aileron chord over local aileron
        Sloc = (c1 + c2) / 2 * (b2 - b1)                    # Surface area of piece
       
        # Determine change in relative airspeed due to yaw
        delta_V = -y_i*m.tan(r*dt)/dt
        V_local = (V_inf + delta_V)*m.cos(beta_nose)
        
        # Determine change in angle of attack due to roll
        roll_induced_alpha = p * y_i / V_local
        alpha_w = alpha_nose + roll_induced_alpha 
        
        # Determine change in angle of attack due to tip vortex
        alpha_i = 0.
        running_alpha_i = True
        while running_alpha_i:
            alpha_e = alpha_w - alpha_i
            Cl, Cd, Cm, xcp = lookup_data(alpha_e, ca_c, da_local,cc)
            alpha_i_new = Cl / (m.pi * AR_w * e_w)
            if (alpha_i_new-alpha_i)/alpha_i_new < 0.01:
                running_alpha_i = False
            alpha_i = alpha_i_new

        alpha_w = alpha_w - alpha_i                 # Angle of Attack as experienced by the piece
        beta_w  = beta_nose                         # Angle of Sideslip as experienced by the piece
        delta_alpha = roll_induced_alpha - alpha_i  # Difference between AoA nose and piece
        
        downwash_angle = 2 * Cl / (m.pi * AR_w)     # Downwash generated by the piece which will be used for the HT
        sidewash_angle = 0.                         # Sidewash generated by the piece which will be used for the VT
        # Lookup Aerodynamic data of the wing based on given imput
        Cl, Cd, Cm, xcp = lookup_data(alpha_w, ca_c, da_local,cc)
        
        # Determine induced drag
        Cdi = Cl * Cl / (m.pi * AR_w * e_w)
        Cd = Cd + Cdi
        
        # Determine moment arm to c.g.
        x_LE = xlemac + (abs(y_i)-y_mac)*m.tan(m.radians(sweep_LE))
        x_i = x_LE + xcp*cc - xcg
        
        # Normal and tangental coefficients (wrt body frame):
        Cn = -Cl*m.cos(alpha_w) - Cd*m.sin(alpha_w)
        Ct = Cl*m.sin(alpha_w) - Cd*m.cos(alpha_w) 
        
        Fn = 0.5 * rho * V_local ** 2 * Sloc * Cn
        Ft = 0.5 * rho * V_local ** 2 * Sloc * Ct

        # Calculate lift and drag in local aerodynamic frame
        Lift = 0.5 * rho * V_local ** 2 * Sloc * Cl
        Drag = 0.5 * rho * V_local ** 2 * Sloc * Cd
        # Calculate lift and drag in general arodynamic frame
        Lift_aero = Lift*m.cos(delta_alpha) - Drag*m.sin(delta_alpha)
        Drag_aero = Drag*m.cos(delta_alpha) + Lift*m.sin(delta_alpha)
        
        # Add everything to the 2D array
        disc_wing_w[i][1] = y_i.magnitude
        disc_wing_w[i][2] = Drag_aero.magnitude
        disc_wing_w[i][3] = Lift_aero.magnitude
        disc_wing_w[i][4] = (Drag * y_i).magnitude
        disc_wing_w[i][5] = (Lift * y_i).magnitude
        disc_wing_w[i][6] = (x_i*Lift).magnitude
        disc_wing_w[i][7] = (x_i*Drag).magnitude
        disc_wing_w[i][8] = downwash_angle
        disc_wing_w[i][9] = sidewash_angle
        
        disc_wing_w[i][11] = Fn.magnitude
        disc_wing_w[i][12] = Ft.magnitude
        disc_wing_w[i][13] = 0.
        disc_wing_w[i][14] = (Fn*y_i).magnitude
        disc_wing_w[i][15] = (Fn*x_i).magnitude
        disc_wing_w[i][16] = (Ft*y_i).magnitude
        disc_wing_w[i][17] = 0.
        disc_wing_w[i][18] = 0.
        disc_wing_w[i][18] = 0.
    print(alpha_w)
    print(disc_wing_w[:,11])    
    # Determine the downwash as function of y-distance to c.g.
    down_wash_func = interpolate.interp1d(disc_wing_w[:,1],disc_wing_w[:,8],'linear')
    side_wash_func = interpolate.interp1d(disc_wing_w[:,1],disc_wing_w[:,9],'linear')

    # Calculate for Horizontal Tail
    for i in range(0, len(khlst)-1):
        de_local = disc_wing_h[i][0]                        # Local aileron deflection of piece
        b1 = khlst[i]                                       # Y boundary left of piece
        b2 = khlst[i+1]                                     # Y boundary right of piece
        y_i = (b1 + b2) / 2                                 # Y centre of piece
        c1 = local_chord(abs(b1), c_r_h, c_t_h, half_b_h)   # Chord left of piece
        c2 = local_chord(abs(b2), c_r_h, c_t_h, half_b_h)   # Chord right of piece
        cc = local_chord(abs(y_i), c_r_h, c_t_h, half_b_h)  # Chord centre of piece
        ce_c = 0.5                                          # percentage aileron chord over local aileron
        Sloc = (c1 + c2) / 2 * (b2 - b1)                    # Surface area of piece
       
        # Determine change in relative airspeed due to yaw
        delta_V = -y_i*m.tan(r*dt)/dt
        V_local = (V_inf + delta_V)*m.cos(beta_nose)
        
        # Determine change in angle of attack due to roll
        roll_induced_alpha = p * y_i / V_local
        downwash = down_wash_func(y_i)
        alpha_h = alpha_nose + roll_induced_alpha - downwash
        
        # Determine change in angle of attack due to tip vortex
        alpha_i = 0.
        running_alpha_i = True
        while running_alpha_i:
            alpha_e = alpha_h - alpha_i
            Cl, Cd, Cm, xcp = lookup_data(alpha_e, ce_c, de_local,cc)
            alpha_i_new = Cl / (m.pi * AR_h * e_h)
            if (alpha_i_new-alpha_i)/alpha_i_new < 0.01:
                running_alpha_i = False
            alpha_i = alpha_i_new

        alpha_h = alpha_h - alpha_i                             # Angle of Attack as experienced by the piece
        beta_h  = beta_nose                                     # Angle of Sideslip as experienced by the piece
        delta_alpha = roll_induced_alpha - alpha_i-downwash     # Difference between AoA nose and piece
        
        # Lookup Aerodynamic data of the wing based on given imput
        Cl, Cd, Cm, xcp = lookup_data(alpha_h, ce_c, de_local,cc)
        
        # Determine induced drag
        Cdi = Cl * Cl / (m.pi * AR_h * e_h)
        Cd = Cd + Cdi
              
        # Normal and tangental coefficients (wrt body frame):
        Cn = -Cl*m.cos(alpha_h) - Cd*m.sin(alpha_h)
        Ct = Cl*m.sin(alpha_h) - Cd*m.cos(alpha_h) 
        
        Fn = 0.5 * rho * V_local ** 2 * Sloc * Cn
        Ft = 0.5 * rho * V_local ** 2 * Sloc * Ct

        # Calculate lift and drag in local aerodynamic frame
        Lift = 0.5 * rho * V_local ** 2 * Sloc * Cl
        Drag = 0.5 * rho * V_local ** 2 * Sloc * Cd
        # Calculate lift and drag in general arodynamic frame
        Lift_aero = Lift*m.cos(delta_alpha) - Drag*m.sin(delta_alpha)
        Drag_aero = Drag*m.cos(delta_alpha) + Lift*m.sin(delta_alpha)
        
        # Add everything to the 2D array
        disc_wing_h[i][1] = Sloc.magnitude
        disc_wing_h[i][2] = Drag_aero.magnitude
        disc_wing_h[i][3] = Lift_aero.magnitude
        disc_wing_h[i][4] = (Drag * y_i).magnitude
        disc_wing_h[i][5] = (Lift * y_i).magnitude
        disc_wing_h[i][6] = (x_i*Lift).magnitude
        disc_wing_h[i][7] = (x_i*Drag).magnitude

        disc_wing_h[i][11] = Fn.magnitude
        disc_wing_h[i][12] = Ft.magnitude
        disc_wing_h[i][13] = 0.
        disc_wing_h[i][14] = (Fn*y_i).magnitude
        disc_wing_h[i][15] = (Fn*l_h).magnitude
        disc_wing_h[i][16] = (Ft*y_i).magnitude
        disc_wing_h[i][17] = 0.
        disc_wing_h[i][18] = 0.
        disc_wing_h[i][19] = 0.



    # Calculate for Vertical Tail
    for i in range(0, len(kvlst)-1):
#        dr_local = disc_wing_v[i][0]
#        b1 = kvlst[i]
#        b2 = kvlst[i+1]
#        z_i = (b1 + b2) / 2
#        c1 = local_chord(abs(b1 - Z_v), c_r_v, c_t_v, b_v)
#        c2 = local_chord(abs(b2 - Z_v), c_r_v, c_t_v, b_v)
#        Sloc = (c1 + c2) / 2 * (b2 - b1)
#        beta = beta_nose + p * z_i / V_inf
#
#        Cl, Cd, Cm, xcp = lookup_data(beta, cr_c, dr_local,1.0*ureg.m)
#
#        Lift = 0.5 * rho * V_inf ** 2 * Sloc * Cl
#        Drag = 0.5 * rho * V_inf ** 2 * Sloc * Cd
#        disc_wing_v[i][1] = Sloc.magnitude
#        disc_wing_v[i][2] = beta.magnitude
#        disc_wing_v[i][3] = Lift.magnitude
#        disc_wing_v[i][4] = c2.magnitude
#        disc_wing_v[i][5] = (Lift.magnitude * z_i.magnitude)
        "---------------------------------------------------------------------"
        dr_local = disc_wing_v[i][0]                        # Local aileron deflection of piece
        b1 = kvlst[i] - Z_v                                 # Y boundary left of piece
        b2 = kvlst[i+1] - Z_v                               # Y boundary right of piece
        z_i = (b1 + b2) / 2                                 # 
        c1 = local_chord(abs(b1), c_r_v, c_t_v, b_v)  # 
        c2 = local_chord(abs(b2), c_r_v, c_t_v, b_v)  # 
        cc = local_chord(abs(z_i), c_r_v, c_t_v, b_v)       # Chord centre of piece
        ca_c = (l_e / cc).magnitude                         # percentage aileron chord over local aileron
        Sloc = (c1 + c2) / 2 * (b2 - b1)                    # Surface area of piece
       
        # Determine change in relative airspeed due to yaw
        V_local = V_inf*m.cos(alpha_nose)
        
        # Determine change in angle of attack due to roll
        roll_induced_beta = p * z_i / V_local
        yaw_induced_beta  = r * l_h / V_local
        sidewash = 0.
        beta_v = round(beta_nose + roll_induced_beta + yaw_induced_beta - sidewash,2)
        
        # Determine change in angle of attack due to tip vortex
        beta_i = 0.
        running_beta_i = True
        while running_beta_i:
            beta_e = beta_v - beta_i
            Cl, Cd, Cm, xcp = lookup_data(beta_e, cr_c, dr_local,cc)
            if Cl == 0.:
                beta_i = 0.
                running_beta_i = False
            else:
                beta_i_new = Cl / (m.pi * AR_v * e_v)
                if (beta_i_new-beta_i)/beta_i_new < 0.01:
                    running_beta_i = False
                beta_i = beta_i_new

        beta_v = beta_v - beta_i                                # Angle of Attack as experienced by the piece
        alpha_v  = alpha_nose                                   # Angle of Sideslip as experienced by the piece
        delta_beta = roll_induced_beta + yaw_induced_beta - beta_i - sidewash      # Difference between AoA nose and piece
        
        # Lookup Aerodynamic data of the wing based on given imput
        Cl, Cd, Cm, xcp = lookup_data(beta_v, ca_c, da_local,cc)
        
        # Determine induced drag
        Cdi = Cl * Cl / (m.pi * AR_v * e_v)
        Cd = Cd + Cdi
        
       
        # Normal and tangental coefficients (wrt body frame):
        Cb = -Cl*m.cos(beta_v) - Cd*m.sin(beta_v)
        Ct = Cl*m.sin(beta_v) - Cd*m.cos(beta_v)

        Ft = 0.5 * rho * V_local ** 2 * Sloc * Ct
        Fb = 0.5 * rho * V_local ** 2 * Sloc * Cb
        
        # Calculate lift and drag in local aerodynamic frame
        Lift = 0.5 * rho * V_local ** 2 * Sloc * Cl
        Drag = 0.5 * rho * V_local ** 2 * Sloc * Cd
        # Calculate lift and drag in general arodynamic frame
        Lift_aero = Lift*m.cos(delta_alpha) - Drag*m.sin(delta_alpha)
        Drag_aero = Drag*m.cos(delta_alpha) + Lift*m.sin(delta_alpha)
        
        # Add everything to the 2D array
        disc_wing_v[i][1] = y_i.magnitude
        disc_wing_v[i][2] = Drag_aero.magnitude
        disc_wing_v[i][3] = Lift_aero.magnitude
        disc_wing_v[i][4] = (Drag * y_i).magnitude
        disc_wing_v[i][5] = (Lift * y_i).magnitude
        disc_wing_v[i][6] = (x_i*Lift).magnitude
        disc_wing_v[i][7] = (x_i*Drag).magnitude
        disc_wing_v[i][8] = downwash_angle
        
        disc_wing_v[i][11] = 0.
        disc_wing_v[i][12] = Ft.magnitude
        disc_wing_v[i][13] = Fb.magnitude
        disc_wing_v[i][14] = 0.
        disc_wing_v[i][15] = 0.
        disc_wing_v[i][16] = 0.
        disc_wing_v[i][17] = (Ft*z_i).magnitude
        disc_wing_v[i][18] = (Fb*z_i).magnitude
        disc_wing_v[i][19] = (Fb*l_v).magnitude



#    # Determine all the Forces and moments caused by the forces    
#    Sum_Dw_y = sum(disc_wing_w[:,4])*ureg.N*ureg.m
#    Sum_Dh_y = sum(disc_wing_h[:,4])*ureg.N*ureg.m
#    Sum_Lw_y = sum(disc_wing_w[:,5])*ureg.N*ureg.m
#    Sum_Lh_y = sum(disc_wing_h[:,5])*ureg.N*ureg.m
#    Sum_Lv_y = sum(disc_wing_v[:,5])*ureg.N*ureg.m
#    Sum_Lw_x = sum(disc_wing_w[:,6])*ureg.N*ureg.m
#    Sum_Dw_x = sum(disc_wing_w[:,7])*ureg.N*ureg.m
#    Sum_Dh_x = 0.*ureg.N*ureg.m
#    Sum_Lw   = sum(disc_wing_w[:,3])*ureg.N
#    Sum_Lh   = sum(disc_wing_h[:,3])*ureg.N
#    Sum_Dw   = sum(disc_wing_w[:,2])*ureg.N
#    Sum_Dh   = sum(disc_wing_h[:,2])*ureg.N
#    Sideforce= sum(disc_wing_v[:,3])*ureg.N
#    
#    # Equations of motion
#    L = - Sum_Lw_y - Sum_Lh_y - Sum_Lv_y 
#    R = (Sum_Dw_y + Sum_Dh_y)*m.cos(beta_nose) -(Sum_Dw_x+Sum_Dh_x)*m.sin(beta_nose) + Sideforce*l_h
#    M = - (Sum_Lh * l_h + Sum_Lw_x)*m.cos(alpha_nose)
#    Fx = T - m.cos(beta_nose)*(Sum_Dw + Sum_Dh + D_fus_gear) - W*m.sin(Theta)
#    Fy = Sideforce*m.cos(beta_nose) - m.sin(beta_nose)*(Sum_Dw + Sum_Dh + D_fus_gear) + W*m.sin(Theta)
#    Fz = - Sum_Lw - Sum_Lh + W*m.cos(Phi)
    
    # New forces and moments
    Sum_Fn     = (sum(disc_wing_w[:,11]) + sum(disc_wing_h[:,11]) + sum(disc_wing_v[:,11])) * Q_("N")
    Sum_Ft     = (sum(disc_wing_w[:,12]) + sum(disc_wing_h[:,12]) + sum(disc_wing_v[:,12])) * Q_("N")
    Sum_Fb     = (sum(disc_wing_w[:,13]) + sum(disc_wing_h[:,13]) + sum(disc_wing_v[:,13])) * Q_("N")
    Sum_Fn_y   = (sum(disc_wing_w[:,14]) + sum(disc_wing_h[:,14]) + sum(disc_wing_v[:,14])) * Q_("N*m")
    Sum_Fn_x   = (sum(disc_wing_w[:,15]) + sum(disc_wing_h[:,15]) + sum(disc_wing_v[:,15])) * Q_("N*m")
    Sum_Ft_y   = (sum(disc_wing_w[:,16]) + sum(disc_wing_h[:,16]) + sum(disc_wing_v[:,16])) * Q_("N*m")
    Sum_Ft_z   = (sum(disc_wing_w[:,17]) + sum(disc_wing_h[:,17]) + sum(disc_wing_v[:,17])) * Q_("N*m")
    Sum_Fb_z   = (sum(disc_wing_v[:,18]) + sum(disc_wing_h[:,18]) + sum(disc_wing_v[:,18])) * Q_("N*m")
    Sum_Fb_y   = (sum(disc_wing_v[:,19]) + sum(disc_wing_h[:,19]) + sum(disc_wing_v[:,19])) * Q_("N*m")    
    #EOM NEW
    # these zeroes are placeholders....
    Mx  = - Sum_Fn_y + Sum_Fb_z
    My  = Sum_Fn_x 
    Mz  = Sum_Ft_y + Sum_Fb_y
    Fx = T + Sum_Ft - W * m.sin(Theta)
    Fy = Sum_Fb + W * m.sin(Phi)
    Fz = Sum_Fn + W * m.cos(Theta) * m.cos(Phi)

    # Kinematic relations
    u_dot = Fx/(mtow)
    u += u_dot*dt
    w_dot = Fz/(mtow)
    w += w_dot*dt
    v_dot = Fy/(mtow)
    v += v_dot*dt
    
    V_inf = np.sqrt(u*u+w*w+v*v)
    
    u_e = V_inf*m.cos(gamma)*m.cos(gamma_lateral)
    w_e = V_inf*m.sin(gamma)
    v_e = V_inf*m.cos(gamma)*m.sin(gamma_lateral)
    
    gamma = m.atan(w_e/u_e)
    gamma_lateral = m.atan(v_e/u_e)
    
    alpha_nose = Theta - gamma 
    beta_nose  = gamma_lateral - Phi
    
    p_dot = Mx / I_xx
    p    += p_dot * dt
    Phi  += p * dt
    
    q_dot  = My / I_yy
    q     += q_dot * dt
    Theta += q * dt
    
    r_dot = Mz / I_zz
    r    += r_dot * dt
    Psi  += r * dt
    
    # update lists for plots
    plst[n]  = p.magnitude
    pdlst[n] = p_dot.magnitude
    qlst[n]  = q.magnitude
    qdlst[n] = q_dot.magnitude
    Fzlst[n] = Fz.magnitude

    n += 1

#    pmax[n_V,1] = max(plst)
#    n_V += 1

# plt.plot(pmax[:,0],pmax[:,1])
# plt.plot(tlst,plst)


print("Finished in:", round(time.time() - t0, 1), "s")



