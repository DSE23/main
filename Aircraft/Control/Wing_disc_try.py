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
t0 = time.time()

# Definitions
def local_chord(z,c_r,c_t,half_b):
    # Calculates the chord at location z(distance from center)
    return c_r-(c_r-c_t)/half_b*z
def lookup_cl(alpha,da):
    dcl_da = da/38.
    return cl_func(math.degrees(alpha)) + dcl_da


# import airfoil lookup tables
f = open("cl_alpha.txt")
cl_alpha = f.readlines()
f.close
cl_alpha_lst = len(cl_alpha)*[0.]
cl_lst = len(cl_alpha)*[0.]
winglst = []

for i in range (0,len(cl_alpha)):
    cl_alpha[i] = cl_alpha[i].strip("\n")
    variables  = cl_alpha[i].split("\t")
    cl_alpha_lst[i] = float(variables[0])
    cl_lst[i]       = float(variables[1])
cl_func = interpolate.interp1d(cl_alpha_lst,cl_lst,'cubic')

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

l_a   = Q_("0.2015 m")                                      # Set aileron length
cabin_width = Geometry.Fuselage.cabin_w + 0.1*ureg.m        # import cabin width 
vt_width    = t_c_v*c_r_v + 0.1*ureg.m                      # import width of VT
rho = IP.rho0                                               # import density

# Setup discritization paramters
n_of_disc_w = 30                # number of parts wing is discretized
n_of_disc_h = 10                # number of parts HT is discretized
n_of_disc_v = 10                # number of parts VT is discretized

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
    kvlst[k]   = 0.312*ureg.m-b_v + bloc_v*k

da    = 30.                     # aileron deflection
alpha_nose = 0.                 # angle of attack of nose
beta_nose  = 0.                 # angle of sideslip of nose
p     = 0./ureg.s               # initial roll rate

Vrange = np.arange(V_s,V_a,1)
pmax = np.transpose(np.vstack([Vrange,np.zeros((1,len(Vrange)))[0]]))

# Init sim
#n_V = 0
#for V in Vrange:
#    V_inf = V*ureg.m/ureg.s
V_inf = 90*ureg.m/ureg.s
t_current = Q_("0 s")
dt = Q_("0.02 s")
t_sim = Q_("0.7 s")
running = True
n = 0
plst = np.zeros((1,int((t_sim/dt).magnitude+1)))[0]
pdlst= np.zeros((1,int((t_sim/dt).magnitude+1)))[0]
tlst = np.arange(0,t_sim.magnitude+dt.magnitude,dt.magnitude)
llst = []
while running:

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
        ca_c = l_a/cc                   # percentage aileron chord over local aileron
        Sloc    = (c1+c2)/2*(b2-b1)     # Surface area of piece
        
        alpha_w   = alpha_nose + p*y_i/V_inf
        alpha_i   = 0.
        for j in range(10):
            alpha_e = alpha_w - alpha_i
            Cl = lookup_cl(alpha_e,da_local)
            alpha_i = Cl/(math.pi*AR_w*e_w)         
        alpha = alpha_w - alpha_i
        downwash_angle = 2*Cl/(math.pi*AR_w)
        
        Lift  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(alpha,da_local)
        Drag  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(alpha,da_local)
        disc_wing_w[i-1][1] = math.degrees(downwash_angle)
        disc_wing_w[i-1][2] = Sloc.magnitude
        disc_wing_w[i-1][3] = Lift.magnitude
        disc_wing_w[i-1][4] = (Drag*y_i).magnitude  # Change to Cd!
        disc_wing_w[i-1][5] = (Lift*y_i).magnitude
    
    for i in range(1,len(khlst)):
        # calculate lift and drag for discretized HT
        b1 = khlst[i-1]
        b2 = khlst[i]
        y_i = (b1+b2)/2
        c1 = local_chord(abs(b1),c_r_h,c_t_h,half_b_h)
        c2 = local_chord(abs(b2),c_r_h,c_t_h,half_b_h)
        Sloc = (c1+c2)/2*(b2-b1)
        alpha = alpha_nose + p*y_i/V_inf
        Lift  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(alpha,0)*1
        Drag  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(alpha,0)*1
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
        beta = beta_nose + p*z_i/V_inf
        Lift  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(beta,0)*1
        Drag  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(beta,0)*1
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
    if t_current.magnitude > t_sim.magnitude:
        running = False
        
    t_current += dt
    n += 1

#    pmax[n_V,1] = max(plst)    
#    n_V += 1

#plt.plot(pmax[:,0],pmax[:,1])
#plt.plot(tlst,plst)
plt.plot(tlst,np.degrees(plst))
plt.plot(tlst,pdlst)

print("Finished in:",round(time.time()-t0,1),"s")