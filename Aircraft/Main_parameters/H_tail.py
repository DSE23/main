# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 21:23:21 2018

@author: jurian
"""

import sys
import math as m
sys.path.append('../')

S_h = 2.629                 # [m^2] Horizontal tail surface
A_h = 2.86                  # Aspect ratio
b_h = m.sqrt(S_h*A_h)       # [m] Span horizontal tail
taper = 0.529               # taper ratio
c_rh = 1.329                # [m] H-tail root chord
c_th = c_rh * taper         # [m] H-tail tip chord
Sweep_h = 0                 # [deg] Sweep H-tail
MAC_h = c_rh*(2/3)*((1+taper+taper**2)/(1+taper))  # [m] Mean aerodynamic chord
S_wet = S_h                 # [m^2] Wetted area
S_e = 1.3145                # Elevator area
c_e = 0.50829               # Elevator chord
delta_e = m.radians(30)     # Max elevator deflection
X_h = 5.27                  # [m] 0.25C location compared to the nose
Z_h = 0.55                  # [m] Distance MAC_h and zero lift line wing

