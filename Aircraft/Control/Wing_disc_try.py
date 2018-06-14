# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7
Author: Jordy van Leeuwen

This code discretizes the main wing, horizontal tail and vertical tail.
Using non-linear EOM and lookup-tables for the aerodynamic properties, the 
manouevre rates and resulting moments are determined. 
"""

import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
import scipy.interpolate as interpolate
import time
from Misc import ureg, Q_     
from Geometry import Geometry
from Aerodynamics import Wing as Awing
from Inertia import Inertia
from Misc import Init_parm as IP
import matplotlib.pyplot as plt
import math
import pandas as pd
t0 = time.time()

# Variables
l_a   = Q_("0.2015 m")                                      # Set aileron length
n_of_disc_w = 30                                            # number of parts wing is discretized
n_of_disc_h = 10                                            # number of parts HT is discretized
n_of_disc_v = 10                                            # number of parts VT is discretized
da    = Q_("30 deg")                                        # aileron deflection
alpha_nose = Q_("0 deg")                                    # angle of attack of nose
beta_nose  = Q_("0 deg")                                    # angle of sideslip of nose
p     = 0./ureg.s                                           # initial roll rate


# Definitions
def local_chord(z,c_r,c_t,half_b):
    # Calculates the chord at location z(distance from center)
    return c_r-(c_r-c_t)/half_b*z


# import airfoil lookup tables
data = pd.read_csv('aerodynamic_data_ms15.dat',' ',header=None).values
def lookup_data(alpha,ca_c,da):
    alpha = math.degrees(alpha)
    da = round(da/5.,0)*5
    ca_c = round(ca_c/0.01,0)*0.01
    
    data_over_da = np.zeros((7,3))
    iterator = 0
    for da_iterate in np.arange(0,31,5):
        row_count_ca = int((ca_c-0.01)*100*51)
        row_count_da = 2550*int(da_iterate/5)
        local_data = data[range(0+row_count_ca+row_count_da,51+row_count_ca+row_count_da),:]
        non_zero_max = max(np.argwhere(local_data[:,0]))[0] # last non-zero row
        local_data = local_data[range(0,non_zero_max+1),:]
        
        if alpha> max(local_data[:,0]):
            if (alpha-max(local_data[:,0]))<1:
                alpha = max(local_data[:,0])
            else:
                raise ValueError("Alpha is out of range")
        if abs(alpha) > max(local_data[:,0]) and alpha<0:
            alpha = -max(local_data[:,0])
        plt.plot(local_data[:,0],local_data[:,1])         
        Cl_local = interpolate.interp1d(local_data[:,0],local_data[:,1],'linear')
        Cd_local = interpolate.interp1d(local_data[:,0],local_data[:,2],'linear')

        if alpha>=0.:
            Cl = Cl_local(alpha)
            Cd = Cd_local(alpha)
        else:
            Cl = -Cl_local(-alpha)
            Cd = Cd_local(-alpha)
        
        data_over_da[iterator,0] = da_iterate
        data_over_da[iterator,1] = Cl
        data_over_da[iterator,2] = Cd    
        iterator +=1
    
    cl_func = interpolate.interp1d(data_over_da[:,0],data_over_da[:,1])
    cd_func = interpolate.interp1d(data_over_da[:,0],data_over_da[:,2])
    
    if da>=0.:
        Cl_final = cl_func(da)*1
        Cd_final = cd_func(da)*1
        #print (data_over_da)
    else:
        Cl_final = -cl_func(-da)*1
        Cd_final = cd_func(-da)*1
    return Cl_final, Cd_final

# Import Aircraft Geometry
I_yy  = Inertia.I_yy
I_xx  = Inertia.I_xx
b_w   = Geometry.Wing.b
S_w   = Geometry.Wing.S
c_r_w = Geometry.Wing.c_r
c_t_w = Geometry.Wing.c_t
AR_w  = Geometry.Wing.A
e_w   = Awing.Oswald_e
b_h   = Geometry.H_tail.b
S_h   = Geometry.H_tail.S
c_r_h = Geometry.H_tail.c_r
c_t_h = Geometry.H_tail.c_t
AR_h  = Geometry.H_tail.A
b_v   = Geometry.V_tail.b
S_v   = Geometry.V_tail.S
c_r_v = Geometry.V_tail.c_r
c_t_v = Geometry.V_tail.c_t
AR_v  = Geometry.V_tail.A
t_c_v = Geometry.V_tail.t_c
Z_v= Geometry.V_tail.Z_v


cabin_width = Geometry.Fuselage.cabin_w + 0.1*ureg.m        # import cabin width 
vt_width    = t_c_v*c_r_v + 0.1*ureg.m                      # import width of VT
rho = IP.rho0                                               # import density


bloc_w = b_w/n_of_disc_w        # span of each station wing
bloc_h = b_h/n_of_disc_h        # Span of each station HT
bloc_v = b_v/n_of_disc_v        # Span of each station VT
half_b_w = b_w/2                # half span wing
half_b_h = b_h/2                # half span HT             

V_s   = IP.V_stall_clean        # stall speed
V_a   = IP.V_a_clean            # manoeuvring speed
V_inf = IP.V_man*ureg.m/ureg.s  # v infinity
  
n_chords_w = int(n_of_disc_w/2)
n_chords_h = int(n_of_disc_h/2)

# Setup lists with station boundaries
kwlst = n_chords_w*[0] + [-cabin_width/2, 0*ureg.m, cabin_width/2] + n_chords_w*[0]
khlst = n_chords_h*[0] + [-vt_width/2, 0*ureg.m, vt_width/2] + n_chords_h*[0]
kvlst = (n_of_disc_v+1)*[0]
for j in range(n_chords_w):
    kwlst[j]    = bloc_w*j - half_b_w
    kwlst[-1-j] = half_b_w - bloc_w*j
for l in range(n_chords_h):
    khlst[l]    = bloc_h*l - half_b_h
    khlst[-1-l] = half_b_h - bloc_h*l
for k in range(n_of_disc_v+1):
<<<<<<< HEAD
    kvlst[k]   = 0.312*ureg.m-b_v + bloc_v*k

da    = 25.                     # aileron deflection
alpha_nose = 0.                 # angle of attack of nose
beta_nose  = 0.                 # angle of sideslip of nose
p     = 0./ureg.s               # initial roll rate
=======
    kvlst[k]   = Z_v-b_v + bloc_v*k
>>>>>>> bf214c8522eb6818eabdc43681adb13a85f3be8f

Vrange = np.arange(V_s,V_a,1)
pmax = np.transpose(np.vstack([Vrange,np.zeros((1,len(Vrange)))[0]]))

# Init sim
#n_V = 0
#for V in Vrange:
#    V_inf = V*ureg.m/ureg.s
V_inf = 90*ureg.m/ureg.s
t_current = Q_("0 s")
dt = Q_("0.01 s")
t_sim = Q_("0.1 s")
running = True
n = 0
plst = np.zeros((1,int((t_sim/dt).magnitude+1)))[0]
pdlst= np.zeros((1,int((t_sim/dt).magnitude+1)))[0]
tlst = np.arange(0,t_sim.magnitude+dt.magnitude,dt.magnitude)
llst = []
while running:
    if (t_current.magnitude+dt.magnitude) > t_sim.magnitude:
        running = False
        
    disc_wing_w = np.zeros((len(kwlst)-1,6))    # 2D array discretized wing
    disc_wing_h = np.zeros((len(khlst)-1,6))    # 2D array discretized HT
    disc_wing_v = np.zeros((len(kvlst)-1,6))    # 2D array discretized VT

    disc_wing_w[(range(n_chords_w)),0] = da                                     # set aileron for postive stations
    disc_wing_w[(range(int(n_of_disc_w+2-n_chords_w),n_of_disc_w+2)),0] = -da   # set aileron for negative stations
    
    for i in range(1,len(kwlst)):
        # calculate lift and drag for discretized wing
        da_local = disc_wing_w[i-1][0]
        b1 = kwlst[i-1]                 # Y boundary left of piece
        b2 = kwlst[i]                   # Y boundary right of piece
        y_i = (b1+b2)/2                 # Y centre of piece
        c1 = local_chord(abs(b1),c_r_w,c_t_w,half_b_w)   # Chord left of piece
        c2 = local_chord(abs(b2),c_r_w,c_t_w,half_b_w)   # Chord right of piece
        cc = local_chord(abs(y_i),c_r_w,c_t_w,half_b_w)
        ca_c = (l_a/cc).magnitude                   # percentage aileron chord over local aileron
        Sloc    = (c1+c2)/2*(b2-b1)     # Surface area of piece
        
        alpha_w   = alpha_nose + p*y_i/V_inf
        alpha_i   = 0.
        for j in range(10):
            alpha_e = alpha_w - alpha_i
            Cl, Cd = lookup_data(alpha_e,ca_c,da_local)
            alpha_i = Cl/(math.pi*AR_w*e_w)         
        alpha = alpha_w - alpha_i
        downwash_angle = 2*Cl/(math.pi*AR_w)
        Cl, Cd = lookup_data(alpha,ca_c,da_local)
        Lift  = 0.5 * rho * V_inf**2 * Sloc * Cl
        Drag  = 0.5 * rho * V_inf**2 * Sloc * Cd
        disc_wing_w[i-1][1] = math.degrees(downwash_angle)
        disc_wing_w[i-1][2] = (p*y_i/V_inf).magnitude
        disc_wing_w[i-1][3] = Lift.magnitude
        disc_wing_w[i-1][4] = (Drag*y_i).magnitude  # Change to Cd!
        disc_wing_w[i-1][5] = (Lift*y_i).magnitude
    print (np.degrees(disc_wing_w[:,2]))
    for i in range(1,len(khlst)):
        # calculate lift and drag for discretized HT
        b1 = khlst[i-1]
        b2 = khlst[i]
        y_i = (b1+b2)/2
        c1 = local_chord(abs(b1),c_r_h,c_t_h,half_b_h)
        c2 = local_chord(abs(b2),c_r_h,c_t_h,half_b_h)
        Sloc = (c1+c2)/2*(b2-b1)
        alpha = alpha_nose #+ p*y_i/V_inf
        Cl, Cd = lookup_data(alpha,0.01,0.)
        Lift  = 0.5 * rho * V_inf**2 * Sloc * Cl
        Drag  = 0.5 * rho * V_inf**2 * Sloc * Cd
        disc_wing_h[i-1][1] = Sloc.magnitude
        disc_wing_h[i-1][2] = y_i.magnitude
        disc_wing_h[i-1][3] = Lift.magnitude
        disc_wing_h[i-1][4] = (Drag*y_i).magnitude  # Change to Cd!
        disc_wing_h[i-1][5] = (Lift.magnitude*y_i.magnitude)   
        
    for i in range(1,len(kvlst)):
        # calculate lift and drag for discretized VT
        b1 = kvlst[i-1]
        b2 = kvlst[i]
        z_i = (b1+b2)/2
        c1 = local_chord(abs(b1-0.312*ureg.m),c_r_v,c_t_v,b_v)
        c2 = local_chord(abs(b2-0.312*ureg.m),c_r_v,c_t_v,b_v)
        Sloc = (c1+c2)/2*(b2-b1)
        beta = beta_nose  + p*z_i/V_inf
        Cl, Cd = lookup_data(beta,0.01,0.)
        Lift  = 0.5 * rho * V_inf**2 * Sloc * Cl
        Drag  = 0.5 * rho * V_inf**2 * Sloc * Cd
        disc_wing_v[i-1][1] = Sloc.magnitude
        disc_wing_v[i-1][2] = beta.magnitude
        disc_wing_v[i-1][3] = Lift.magnitude
        disc_wing_v[i-1][4] = c2.magnitude  # Change to Cd!
        disc_wing_v[i-1][5] = (Lift.magnitude*z_i.magnitude) 
        
    L = (- sum(disc_wing_w[:,5]) -sum(disc_wing_h[:,5]) - sum(disc_wing_v[:,5])) * ureg.N*ureg.m
    R = (sum(disc_wing_w[:,4]) + sum(disc_wing_h[:,4])) * ureg.N*ureg.m
    
    p_dot = L/I_xx    
    p += p_dot*dt
    
    # update lists for plots
    plst[n] = p.magnitude
    pdlst[n]= p_dot.magnitude
    llst.append(L)
    
    # terminate simulation at specified end time

        
    t_current += dt
    n += 1

#    pmax[n_V,1] = max(plst)    
#    n_V += 1

#plt.plot(pmax[:,0],pmax[:,1])
#plt.plot(tlst,plst)
plt.plot(tlst,np.degrees(plst))
plt.plot(tlst,pdlst)

plt.show()

print("Finished in:",round(time.time()-t0,1),"s")