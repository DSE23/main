# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 20:56:57 2018

@author: jurian
"""

import sys
import math as m
sys.path.append('../')

S = 11.74                   # [m^2] Wing Surface
A = 5.5                     # Aspect Ratio
b = m.sqrt(S*A)             # [m] Wing Span
taper = 0.45                # Taper ratio
c_r = 2.015                 # Root chord
c_t = c_r * taper           # Tip chord
Sweep_25 = 0                # [deg] Quarter chord sweep
Sweep_50 = m.degrees(m.atan(m.tan(m.radians(Sweep_25))-(4/A) *
                     ((0.5-0.25)*(1-taper)/(1+taper))))   # [deg] Half chord sweep
Sweep_75 = m.degrees(m.atan(m.tan(m.radians(Sweep_25))-(4/A) *
                     ((0.75-0.25)*(1-taper)/(1+taper))))   # [deg] 3/4 chord sweep
Sweep_LE = m.degrees(m.atan(m.tan(m.radians(Sweep_25))-(4/A) *
                     ((0.0-0.25)*(1-taper)/(1+taper))))   # [deg] LE chord sweep
Dihedral = 0.0              # [deg] Dihedral angle
MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))   # [m] Mean aerodynamic chord
S_wet = S                   # Wetted wing area
S_a = 2.677                 # Aileron area
c_a = 0.3828                # Aileron chord
delta_a = m.radians(30)     # [rad] max aileron deflection
delta_CL_max_a = 0.8267     # Max lift coeff difference due to aileron deflection
