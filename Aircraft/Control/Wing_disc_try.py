# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:19:31 2018

@author: jurian
"""
import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path
import numpy as np
import scipy.interpolate as interpolate
import time
from Misc import ureg, Q_     
from Geometry import Geometry
from Inertia import Inertia
from Misc import Init_parm as IP
import matplotlib.pyplot as plt
import math
t0 = time.time()


def local_chord(z,c_r,c_t,half_b):
    # Calculates the chord at location z(distance from center)
    return c_r-(c_r-c_t)/half_b*z


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
def lookup_cl(alpha,da):
    dcl_da = da/38.
    return cl_func(math.degrees(alpha)) + dcl_da


I_yy  = Inertia.I_yy                         # import MOI around y
I_xx  = Inertia.I_xx                         # import MOI around x
b_w   = Geometry.Wing.b                      # import Span wing
b_h   = Geometry.H_tail.b                    # import Span HT
S_w   = Geometry.Wing.S                      # import Surface area
S_h   = Geometry.H_tail.S
c_r_w = Geometry.Wing.c_r                    # import Root Chord
c_t_w = Geometry.Wing.c_t                    # calculate Tip Chord
c_r_h = Geometry.H_tail.c_r
c_t_h = Geometry.H_tail.c_t

t_c_v = Geometry.V_tail.t_c
c_v   = Geometry.V_tail.c_r


cabin_width = Geometry.Fuselage.cabin_w + 0.1*ureg.m         # import cabin width 
vt_width    = t_c_v*c_v + 0.1*ureg.m
rho = IP.rho0                               # import density

n_of_disc_w = 40                # number of parts wing is discretized
n_of_disc_h = 10                # number of parts HT is discretized
bloc_w = b_w/n_of_disc_w        # span of each station wing
bloc_h = b_h/n_of_disc_h        # Span of each station HT
half_b_w = b_w/2                # half span wing
half_b_h = b_h/2                # half span HT
V_s   = IP.V_stall_clean        # stall speed
V_a   = IP.V_a_clean            # manoeuvring speed
V_inf = IP.V_man*ureg.m/ureg.s                # v infinity
  

n_chords_w = int(n_of_disc_w/2)
n_chords_h = int(n_of_disc_h/2)

kwlst = n_chords_w*[0] + [-cabin_width/2, 0*ureg.m, cabin_width/2] + n_chords_w*[0]      # list containing all station boundaries
khlst = n_chords_h*[0] + [-vt_width/2, 0*ureg.m, vt_width/2] + n_chords_h*[0]
for j in range(n_chords_w):
    kwlst[j]    = bloc_w*j - half_b_w
    kwlst[-1-j] = half_b_w - bloc_w*j
for l in range(n_chords_h):
    khlst[l]    = bloc_h*l - half_b_h
    khlst[-1-l] = half_b_h - bloc_h*l
disc_wing_w = np.zeros((len(kwlst)-1,6))           # 2D array discretized wing: Sloc - b_centre - Delta lift -  Delta Lift * b_centre
disc_wing_h = np.zeros((len(khlst)-1,6))

da    = 30                      # aileron deflection
alpha_nose = 0.                      # angle of attack
beta  = 0.                      # angle of sideslip
p     = 0./ureg.s                      # roll

Vrange = np.arange(V_s,V_a,1)
pmax = np.transpose(np.vstack([Vrange,np.zeros((1,len(Vrange)))[0]]))

# Init sim
n_V = 0
for V in Vrange:
    t_current = Q_("0 s")
    dt = Q_("0.05 s")
    t_sim = Q_("1.0 s")
    
    plst = np.zeros((1,int((t_sim/dt).magnitude+1)))[0]
    tlst = np.arange(0,t_sim.magnitude+dt.magnitude,dt.magnitude)
    
    V_inf = V*ureg.m/ureg.s
    running = True
    n = 0
    while running:
        disc_wing_w[(range(n_chords_w)),0] = da                                # set aileron for postive stations
        disc_wing_w[(range(int(n_of_disc_w+2-n_chords_w),n_of_disc_w+2)),0] = -da  # set aileron for negative stations
        
        for i in range(1,len(kwlst)):
            da_local = disc_wing_w[i-1][0]
            b1 = kwlst[i-1]              # Y boundary left of piece
            b2 = kwlst[i]                # Y boundary right of piece
            y_i = (b1+b2)/2              # Y centre of piece
            c1 = local_chord(abs(b1),c_r_w,c_t_w,half_b_w)   # Chord left of piece
            c2 = local_chord(abs(b2),c_r_w,c_t_w,half_b_w)   # Chord right of piece
            Sloc  = (c1+c2)/2*(b2-b1)    # Surface area of piece
            alpha =   alpha_nose + p*y_i/V_inf
            Lift  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(alpha,da_local)
            Drag  = 0.5 * rho * V_inf**2 * Sloc * lookup_cl(alpha,da_local)
            disc_wing_w[i-1][1] = Sloc.magnitude
            disc_wing_w[i-1][2] = y_i.magnitude
            disc_wing_w[i-1][3] = Lift.magnitude
            disc_wing_w[i-1][4] = Drag.magnitude  # Change to Cd!
            disc_wing_w[i-1][5] = (Lift*y_i).magnitude
    
        
        for i in range(1,len(khlst)):
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
            disc_wing_h[i-1][4] = Drag.magnitude  # Change to Cd!
            disc_wing_h[i-1][5] = (Lift.magnitude*y_i.magnitude)        
            
        L = - sum(disc_wing_w[:,5])*ureg.N*ureg.m -sum(disc_wing_h[:,5])*ureg.N*ureg.m
        p_dot = L/I_xx    
        p += p_dot*dt
        
        plst[n] = (p.magnitude)
    
        if t_current.magnitude > t_sim.magnitude:
            running = False
            
        t_current += dt
        n += 1

    pmax[n_V,1] = max(plst)    
    n_V += 1
print("Finished in:",round(time.time()-t0,1),"s")
plt.plot(pmax[:,0],pmax[:,1])
