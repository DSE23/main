# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:19:31 2018

@author: jurian
"""

import sys
sys.path.append('../')    # This makes sure the parent directory gets added to the system path

import numpy as np
from Misc import ureg, Q_      # Imports the unit registry fron the Misc folder
from Geometry import Wing
from Misc import Init_parm as IP


b = Wing.s                      # import Span
c_r = Wing.ChordR               # import Root Chord
taper = Wing.t                  # import Taper
c_t = c_r*taper                 # calculate Tip Chord


def local_chord(z):  # Calculates the chord at location z(distance from center)
    chord = c_r-(c_r-c_t)/half_b*z
    return chord


n_of_disc = 20                  # Number of parts wing is discretized
bloc = b/n_of_disc              # Span of each piece
half_b = b/2
V_inf = IP.V_man

Slst = n_of_disc * [0]
b_centre_lst = n_of_disc * [0]

for i in range(n_of_disc):
    b1 = bloc*i - half_b        # Z boundary left
    b2 = bloc*(i+1) - half_b    # Z boundary right
    b_centre = (b1+b2)/2        # Z centre of piece
    c1 = local_chord(abs(b1))   # Left chord
    c2 = local_chord(abs(b2))   # Right chord
    Sloc = ((c1+c2)/2)*(b2-b1)  # Surface area of piece
    Slst[i] = Sloc
    b_centre_lst[i] = b_centre
    
print (sum(Slst))