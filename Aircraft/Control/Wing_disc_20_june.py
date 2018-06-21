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
stdout_old = sys.stdout
sys.stdout = open(os.devnull, 'w')
sys.path.append('../')              # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
from Geometry import Geometry
from Aerodynamics import Wing as Awing
from Inertia import Inertia
from Performance import Performance
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import pandas as pd
import math as m
import time
sys.stdout = stdout_old
t0 = time.time()

np.seterr(all='raise')

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
alpha_nose = Q_("0. rad") # angle of attack of nose
beta_nose  = Q_("0. rad")   # angle of sideslip of nose
V_inf = Q_("60 m/s")     # V infinity
t_current = Q_("0.0 s")       # Start time of sim
dt = Q_("0.05 s")           # Time step of sim
t_end = Q_("10. s")         # End time of sim
l_h = Q_("3.6444 m")        # Tail arm ac-ac horizontal
l_v = Q_("3.7 m")           # Tail arm ac-ac vertical
p = Q_("0. 1/s")            # initial roll rate  [rad/s]
q = Q_("0. 1/s")            # initial pitch rate [rad/s]
r = Q_("0. 1/s")            # initial yaw rate   [rad/s]
Phi   = Q_("0. rad")        # Initial euler angle around x-axis
Psi   = Q_("0. rad")        # Initial euler angle around z-axis
Theta = Q_("0. rad")        # Initial euler angle around y-axis
lin_ran_alpha = 4          # Linear range of angle of attack and elevator defl.
w = Q_("0. m/s")
u = V_inf
v = Q_("0. m/s")
gamma = Q_("0. rad")
Xi = Q_("0 rad")
mtow = Geometry.Masses.W_MTOW
g0 = Performance.g0.magnitude * Q_("m/s**2")

# Import Aircraft Geometry
I_yy = Inertia.I_yy
I_xx = Inertia.I_xx
I_zz = Inertia.I_zz
I_xz = Inertia.I_xz
I_star = I_xx*I_zz-I_xz**2

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
MAC_htail = Geometry.H_tail.MAC
X_h = Geometry.CG.X_htail
X_w = Geometry.CG.X_wing

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
def trimming(u,ca_c, ce_c):
    # Calculates trim condition
    Cl1, Cdhde, Cmhde, Xcphde =  (lookup_data(0,ce_c,0))                   # Change to H_tail
    Cl2, Cd2hde, Cm2hde, Xcp2hde =  (lookup_data(0,ce_c,lin_ran_alpha))       #       ""    
    dCl_de = (Cl2-Cl1)/(lin_ran_alpha)                          # Elevator effectiveness
    dCd_de = (Cd2hde-Cdhde)/lin_ran_alpha
    dCm_de = (Cm2hde-Cmhde)/lin_ran_alpha
    Cla1w, Cda1w, Cma1w, Xcpa1w = (lookup_data(0., ca_c, 0))                 
    Cla2w, Cda2w, Cma2w, Xcpa2w = (lookup_data(m.radians(lin_ran_alpha), ca_c, 0))
    dCl_alpha_w = (Cla2w-Cla1w)/m.radians(lin_ran_alpha)     # Cl alpha calculation of Wing
    dCd_alpha_w = (Cda2w-Cda1w)/m.radians(lin_ran_alpha)     # Cd alpha calc of wing
    dCm_alpha_w = (Cma2w - Cma1w)/m.radians(lin_ran_alpha)   # Cm alpha calc of wing
    Cla1h, Cda1h, Cma1h, Xcpa1h = (lookup_data(0., ce_c, 0))                 
    Cla2h, Cda2h, Cma2h, Xcpa2h = (lookup_data(m.radians(lin_ran_alpha), ce_c, 0))
    dCl_alpha_h = (Cla2h - Cla1h)/m.radians(lin_ran_alpha)     # Cl alpha calc. of H_tail
    dCd_alpha_h = (Cda2h - Cda1h)/m.radians(lin_ran_alpha)     # Cd alpha calc. of H_tail
    dCm_alpha_h = (Cma2h - Cma1h)/m.radians(lin_ran_alpha)     # Cm alpha calc. of H_tail
    q_Sw = S_w * 0.5 * rho * (u**2)          # Dynamic press. times Wing surface
    q_Sh = S_h * 0.5 * rho * (u**2)          # Dyn. press. times H_tail surface
    l_w = X_w - xcg                                 # Wing arm
    l_h = X_h - xcg                                 # Tail arm
    trim_mat = np.matrix([[(q_Sw*(Cda1w-dCl_alpha_w)).magnitude, (q_Sh*(Cdhde-dCl_alpha_h)).magnitude,
                         (-q_Sh*dCl_de).magnitude],
                        [(q_Sw*((Cda1w-dCl_alpha_w)*l_w+dCm_alpha_w*MAC)).magnitude, (q_Sh*((Cda1h-dCl_alpha_h)*l_h+dCm_alpha_h*MAC_htail)).magnitude,
                         (q_Sh*(dCm_de*MAC_htail - dCl_de*l_h)).magnitude],
                         [1-2*dCl_alpha_w/(m.pi*AR_w * e_w), -1, 0]])
    trim_mat2 = np.matrix([[((-mtow*g0).magnitude*np.cos(Theta))],
                          [0],
                          [0]])
    trim_cond = np.linalg.solve(trim_mat, trim_mat2)
    alpha_t = trim_cond[0] + (dCl_alpha_w * trim_cond[0])/(m.pi * AR_w * e_w)
#    print((trim_cond[0]))
#    print((trim_cond[1]))
#    print((trim_cond[0]-2*trim_cond[0]*dCl_alpha_w/(m.pi*AR_w * e_w)))
    
    alpha_t = m.degrees(alpha_t)
    de_t = trim_cond[2]
    return alpha_t, de_t


# import airfoil lookup tables
data = pd.read_csv('aerodynamic_data_ms15.dat', ' ', header=None).values

def lookup_data(alpha, ca_c, da):
    # Looksup and interpolates Cl and Cd based on alpha, ca_c and da
    alpha = round(alpha,5)
    alpha = m.degrees(alpha)
    indexca_c = int(100*ca_c)-1
    if da == 0.:
        indexda = abs(da)//5
        localdata = data[int(indexda*50*51+indexca_c*51):int(indexda*50*51+(indexca_c+1)*51),:]
        non_zero_max = max(np.argwhere(localdata[:, 0]))[0]  # last non-zero row
        localdata = localdata[:non_zero_max+1,:]
    else:
        index1da = abs(da)//5
        index2da = abs(da)//5 + 1
        index3da = 0
        localdata1 = data[int(index1da*50*51+indexca_c*51):int(index1da*50*51+(indexca_c+1)*51),:]
        localdata2 = data[int(index2da*50*51+indexca_c*51):int(index2da*50*51+(indexca_c+1)*51),:]
        localdata3 = data[int(index3da*50*51+indexca_c*51):int(index3da*50*51+(indexca_c+1)*51),:]
        non_zero_max1 = max(np.argwhere(localdata1[:, 0]))[0]  # last non-zero row
        non_zero_max2 = max(np.argwhere(localdata2[:, 0]))[0]  # last non-zero row
        non_zero_max3 = max(np.argwhere(localdata3[:, 0]))[0]  # last non-zero row
        non_zero_max = min((non_zero_max1,non_zero_max2,non_zero_max3))
        
        localdata1 = localdata1[:non_zero_max+1,:]
        localdata2 = localdata2[:non_zero_max+1,:]    
        localdata3 = localdata3[:non_zero_max+1,:]   
        
        localdata_alpha = localdata2[:,0]
        localdata_cl    = (localdata2[:,1] - localdata1[:,1])/5 *abs(da) + localdata3[:,1]
        localdata_cd    = (localdata2[:,2] - localdata1[:,2])/5 *abs(da) + localdata3[:,2]
        localdata_cm    = (localdata2[:,4] - localdata1[:,4])/5 *abs(da) + localdata3[:,4]
    
        localdata = np.array([localdata_alpha,
                              localdata_cl,
                              localdata_cd,
                              localdata_cd,
                              localdata_cm])
        localdata = np.transpose(localdata)
    
    Cl_local = interpolate.interp1d(localdata[:,0], localdata[:,1], 'linear', fill_value='extrapolate')
    Cd_local = interpolate.interp1d(localdata[:,0], localdata[:,2], 'linear', fill_value='extrapolate')
    Cm_local = interpolate.interp1d(localdata[:,0], localdata[:,4], 'linear', fill_value='extrapolate')
    if da >= 0:
        Cl = Cl_local(alpha)
        Cd = Cd_local(alpha)
        Cm = Cm_local(alpha)
    else:
        Cl = -Cl_local(-alpha)
        Cd = Cd_local(-alpha)
        Cm = -Cm_local(-alpha)
    
    Cn_lookup = m.cos(m.radians(alpha))*Cl + m.sin(m.radians(alpha))*Cd
    if Cn_lookup !=0:
        xcp = 0.25-Cm/Cn_lookup
    else:
        xcp = 0.25
    return Cl, Cd, Cm, xcp
    
# Import other forces
T = Q_("1000 N")
T_g_f = 0. # UPDATE REQUIRED
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
running = True
n = 0
# empty arrays for plotting in the end:
plst  = np.zeros((1, int((t_end / dt).magnitude)))[0]
alst  = np.zeros((1, int((t_end / dt).magnitude)))[0]
pdlst = np.zeros((1, int((t_end / dt).magnitude)))[0]
qlst  = np.zeros((1, int((t_end / dt).magnitude)))[0]
qdlst = np.zeros((1, int((t_end / dt).magnitude)))[0]
Fzlst = np.zeros((1, int((t_end / dt).magnitude)))[0]
tlst  = np.arange(0, t_end.magnitude, dt.magnitude)
t2lst = np.zeros((1, int((t_end / dt).magnitude)))[0]
alpha_nose,de = trimming(V_inf,0.1,0.5)
alpha_nose = m.radians(alpha_nose)
Theta = alpha_nose
de = de[0,0]


#raise ValueError("breakie breakie")
for t_current in np.arange(0,(t_end).magnitude,dt.magnitude):
    t_start_loop = time.time()
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
            Cl, Cd, Cm, xcp = lookup_data(alpha_e, ca_c, da_local)
            if abs(Cl) <0.01:
                alpha_i_new = 0.
                running_alpha_i = False
            else:
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
        Cl, Cd, Cm, xcp = lookup_data(alpha_w, ca_c, da_local)
        
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
            Cl, Cd, Cm, xcp = lookup_data(alpha_e, ce_c, da_local)
            if abs(Cl) <0.01:
                alpha_i_new = 0.
                running_alpha_i = False
            else:
                alpha_i_new = Cl / (m.pi * AR_h * e_h)
                if (alpha_i_new-alpha_i)/alpha_i_new < 0.01:
                    running_alpha_i = False
                alpha_i = alpha_i_new

        alpha_h = alpha_h - alpha_i                             # Angle of Attack as experienced by the piece
        beta_h  = beta_nose                                     # Angle of Sideslip as experienced by the piece
        delta_alpha = roll_induced_alpha - alpha_i-downwash     # Difference between AoA nose and piece
        

        # Lookup Aerodynamic data of the wing based on given imput
        Cl, Cd, Cm, xcp = lookup_data(alpha_h, ce_c, de_local)
        
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
            Cl, Cd, Cm, xcp = lookup_data(beta_e, cr_c, dr_local)
            if abs(Cl) <0.01:
                alpha_i_new = 0.
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
        Cl, Cd, Cm, xcp = lookup_data(beta_v, ca_c, da_local)
        
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
    Fx = T + Sum_Ft + T_g_f - W * m.sin(Theta)
    Fy = Sum_Fb + W * m.sin(Phi) * m.cos(Theta)
    Fz = Sum_Fn + W * m.cos(Theta) * m.cos(Phi)

    Mx = Sum_Fn_y + Sum_Fb_z
    My = Sum_Fn_x
    Mz = -Sum_Ft_y + Sum_Fb_y


    # Kinematic relations
    u_dot = Fx/(mtow) - q*w + r*v
    v_dot = Fy/(mtow) - r*u + p*w
    w_dot = Fz/(mtow) - p*v + q*u
    
    u += u_dot*dt
    w += w_dot*dt
    v += v_dot*dt
    
    V_inf = np.sqrt(u*u+w*w+v*v)
    
    u_e = u*(m.cos(Theta)*m.cos(Psi)) +\
          v*(m.sin(Phi)*m.sin(Theta)*m.cos(Psi) - m.cos(Phi)*m.sin(Psi)) +\
          w*(m.cos(Phi)*m.sin(Theta)*m.cos(Psi) + m.sin(Phi)*m.sin(Psi))
    v_e = u*(m.cos(Theta)*m.sin(Psi)) +\
          v*(m.sin(Phi)*m.sin(Theta)*m.sin(Psi) + m.cos(Phi)*m.cos(Psi)) +\
          w*(m.cos(Phi)*m.sin(Theta)*m.sin(Psi))
    w_e = u*(m.sin(Theta)) +\
          v*(m.sin(Phi)*m.cos(Theta)) +\
          w*(m.cos(Phi)*m.cos(Theta))
    
    gamma = m.atan(-w_e/(np.sqrt(u_e**2+v_e**2)))
    Xi = m.atan(v_e/u_e)
    
    p_dot = I_zz/I_star * Mx + I_xz/I_star * Mz +\
            ((I_xx-I_yy+I_zz)*I_xz)/I_star * p * q +\
            ((I_yy-I_zz)*I_zz-I_xz**2)/I_star * q * r
    q_dot = My/I_yy + I_xz/I_yy * (r**2 - p**2) +\
            (I_zz-I_xx)/I_yy * p * r
    r_dot = I_xz/I_star * Mx + I_xx/I_star * Mz +\
            ((I_xx-I_yy)*I_xx + I_xz**2)/I_star * p * q +\
            ((-I_xx+I_yy-I_zz)*I_xz)/I_star * p *r
    
    p += p_dot * dt
    q += q_dot * dt
    r += r_dot * dt
    
    Phi += (p + q*m.sin(Phi)*m.tan(Theta)  + r*m.cos(Phi)*m.tan(Theta)) * dt
    Theta += (q*m.cos(Phi) - r*m.sin(Phi)) * dt
    Psi += (q*m.sin(Phi)/m.cos(Theta) + r*m.cos(Phi)/m.cos(Theta)) * dt
    
    
    alpha_nose= Theta - gamma
    beta_nose = Psi - Xi
    
    # update lists for plots
    plst[n]  = p.magnitude
    pdlst[n] = p_dot.magnitude
    qlst[n]  = q.magnitude
    qdlst[n] = q_dot.magnitude
    Fzlst[n] = Fz.magnitude
    alst[n] = alpha_nose.magnitude

    
    t_end_loop = time.time()
    t2lst[n] = t_end_loop-t_start_loop
    n += 1
    #print("Lw:",sum(disc_wing_w[:,11]))
plt.plot(tlst,np.degrees(alst),label="alpha")
plt.plot(tlst,np.degrees(qlst),label="pitch")
plt.legend()
plt.show()
# plt.plot(pmax[:,0],pmax[:,1])
# plt.plot(tlst,plst)
#plt.plot(tlst,Fzlst)

print("Average Loop Time:",round(np.average(t2lst),2),"s")
print("Finished in:", round(time.time() - t0, 1), "s")
