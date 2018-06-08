# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:19:31 2018

@author: jurian
"""

import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path

import numpy as np
import scipy as sp
from Misc import ureg, Q_      # Imports the unit registry fron the Misc folder
from Geometry import Wing
from Geometry import Fuselage
from Inertia import Inertia
from Misc import Init_parm as IP

f = open("cl_alpha.txt")
cl_alpha = f.readlines()
f.close

g = open("cd_alpha.txt")
cd_alpha = g.readlines()
g.close
cl_alpha_lst = len(cl_alpha)*[0.]
cl_lst = len(cl_alpha)*[0.]


for i in range (0,len(cl_alpha)):
    cl_alpha[i] = cl_alpha[i].strip("\n")
    variables  = cl_alpha[i].split("\t")
    cl_alpha_lst[i] = float(variables[0])
    cl_lst[i]       = float(variables[1])
cl_func = sp.interpolate.interp1d(cl_alpha_lst,cl_lst,'cubic')


Iyy = Inertia.Iyy_aircraft        
b = Wing.s                      # import Span
c_r = Wing.ChordR               # import Root Chord
taper = Wing.t                  # import Taper
c_t = c_r*taper                 # calculate Tip Chord
cabin_w = Fuselage.cabin_w      # 0.6 cabin width + 2*0.1 structure 

def local_chord(z):  # Calculates the chord at location z(distance from center)
    chord = c_r-(c_r-c_t)/half_b*z
    return chord

n_of_disc = 50                  # Number of parts wing is discretized
bloc = b/n_of_disc              # Span of each piece
half_b = b/2
V_inf = IP.V_man
  

n_chords = int(n_of_disc/2)
klst = n_chords*[0] + [-cabin_w/2, 0*ureg.m, cabin_w/2] + n_chords*[0]
for j in range(n_chords):
    klst[j]=bloc*j - half_b
    klst[-1-j] = half_b - bloc*j

Slst = (len(klst)-1)*[0]
b_centrelst = (len(klst)-1)*[0]
disc_wing = np.zeros((len(klst)-1,3))

for i in range(1,len(klst)):
    b1 = klst[i-1]              # Z boundary left of piece
    b2 = klst[i]                # Z boundary right of piece
    b_centre = (b1+b2)/2        # Z centre of piece
    c1 = local_chord(abs(b1))   # Chord left of piece
    c2 = local_chord(abs(b2))   # Chord right of piece
    Sloc = (c1+c2)/2*(b2-b1)    # Surface area of piece
    disc_wing[i-1][0] = Sloc.magnitude
    disc_wing[i-1][1] = b_centre.magnitude
    

    
