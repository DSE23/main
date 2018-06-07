# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 21:40:11 2018

@author: jurian
"""

import sys
import math as m
sys.path.append('../')

S = 1.08                                    # [m^2] V-tail surface area
A = 0.99                                    # Aspect ratio
b = m.sqrt(S*A)                             # [m] Span V-tail
taper = 0.33                                # Taper ratio
c_r = 1.660                                 # [m] V-tail root chord
c_t = c_r * taper                           # [m] V-tail tip chord
Sweep = 0                                   # [deg] Sweep V-tail
MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))  # [m] Mean aerodynamic chord
S_wet = S                                   # [m^2] Wetted area
S_r = 0.5726                                # Rudder area
c_r = 0.553                                 # Rudder chord
delta_r = m.radians(30)                     # Max rudder deflection
X_v = 5.70                                  # [m] 0.25C location compared to the nose
Z_v = 0.52                                  # [m] Distance MAC_h and zero lift line wing
