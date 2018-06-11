# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:19:31 2018

@author: jurian
"""

import sys
import numpy as np
import scipy as sp
import time
from Misc import ureg, Q_      # Imports the unit registry fron the Misc folder
from Geometry import Geometry
from Inertia import Inertia
from Misc import Init_parm as IP
t0 = time.time()

sys.path.append('../')    # This makes sure the parent directory gets added to the system path
def local_chord(z):  # Calculates the chord at location z(distance from center)
    chord = c_r-(c_r-c_t)/half_b*z
    return chord

f = open("cl_alpha.txt")
cl_alpha = f.readlines()
f.close
cl_alpha_lst = len(cl_alpha)*[0.]
cl_lst = len(cl_alpha)*[0.]


for i in range (0,len(cl_alpha)):
    cl_alpha[i] = cl_alpha[i].strip("\n")
    variables  = cl_alpha[i].split("\t")
    cl_alpha_lst[i] = float(variables[0])
    cl_lst[i]       = float(variables[1])
cl_func = sp.interpolate.interp1d(cl_alpha_lst,cl_lst,'cubic')


I_yy = Inertia.I_yy        
I_xx = Inertia.I_xx
b   = Geometry.Wing.b                      # import Span
S   = Geometry.Wing.S
c_r = Geometry.Wing.c_r               # import Root Chord
taper   = Geometry.Wing.taper                  # import Taper
c_t = c_r*taper                 # calculate Tip Chord
cabin_w = Geometry.Fuselage.cabin_w      # 0.6 cabin width + 2*0.1 structure 



n_of_disc = 40                  # Number of parts wing is discretized
bloc = b/n_of_disc              # Span of each piece
half_b = b/2                    # Half span
V_s   = IP.V_stall_clean        # Stall speed
V_a   = IP.V_a_clean            # Manoeuvring speed
V_inf = IP.V_man                # V infinity
  

n_chords = int(n_of_disc/2)
klst = n_chords*[0] + [-cabin_w/2, 0*ureg.m, cabin_w/2] + n_chords*[0]      # list containing all station boundaries
for j in range(n_chords):
    klst[j]=bloc*j - half_b
    klst[-1-j] = half_b - bloc*j

disc_wing = np.zeros((len(klst)-1,5))           # 2D array discretized wing: Sloc - b_centre - Delta lift -  Delta Lift * b_centre
da = 30

rho = IP.rho0
alpha = 0.
p = 0.

t_current = 0
dt = 0.05
running = True
while running:
    disc_wing[(range(n_chords)),0] += da                                # set aileron for postive stations
    disc_wing[(range(int(n_of_disc+2-n_chords),n_of_disc+2)),0] += -da  # set aileron for negative stations
    for i in range(1,len(klst)):
        b1 = klst[i-1]              # Z boundary left of piece
        b2 = klst[i]                # Z boundary right of piece
        b_centre = (b1+b2)/2        # Z centre of piece
        c1 = local_chord(abs(b1))   # Chord left of piece
        c2 = local_chord(abs(b2))   # Chord right of piece
        Sloc = (c1+c2)/2*(b2-b1)    # Surface area of piece
        disc_wing[i-1][1] = Sloc.magnitude
        disc_wing[i-1][2] = b_centre.magnitude
        disc_wing[i-1][3] = (0.5 * rho * V_inf**2 * Sloc * cl_func(alpha)*1).magnitude
        disc_wing[i-1][4] = (disc_wing[i-1][2]*b_centre).magnitude
        
    L = sum(disc_wing[:,4])*ureg.N*ureg.m
    p_dot = L/I_xx    
    p += p_dot
    
    
    if t_current >1.:
        running = False
        
    t_current += dt

    
print("Finished in:",round(time.time()-t0,1),"s")