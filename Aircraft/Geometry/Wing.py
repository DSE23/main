# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 20:56:57 2018

@author: jurian
Changed by boris 
DONT FORGET TO ADD THE UNITS TO YOUR VARIABLES......
"""

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

import numpy as np
import scipy as sp
from scipy import interpolate
import math as m
from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder


S = Q_('11.74 m**2')                   # [m^2] Wing Surface
A = 5.5                     # Aspect Ratio
b = np.sqrt(S*A)             # [m] Wing Span
taper = 0.45                # Taper ratio
c_r = Q_('2.015 m')                 # Root chord
c_t = c_r * taper           # Tip chord
Sweep_25 = 0                # [deg] Quarter chord sweep
Sweep_25 *= Q_('deg')
Sweep_50 = m.degrees(m.atan(m.tan(m.radians(Sweep_25))-(4/A) *
                     ((0.5-0.25)*(1-taper)/(1+taper))))   # [deg] Half chord sweep
Sweep_50 *= Q_('deg')
Sweep_75 = m.degrees(m.atan(m.tan(m.radians(Sweep_25))-(4/A) *
                     ((0.75-0.25)*(1-taper)/(1+taper))))   # [deg] 3/4 chord sweep
Sweep_75 *= Q_('deg')
Sweep_LE = m.degrees(m.atan(m.tan(m.radians(Sweep_25))-(4/A) *
                     ((0.0-0.25)*(1-taper)/(1+taper))))   # [deg] LE chord sweep
Sweep_LE *= Q_('deg')
Dihedral = Q_('0.0 deg')             # [deg] Dihedral angle
MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))   # [m] Mean aerodynamic chord
S_wet = S                   # Wetted wing area
S_a = Q_('2.677 m**2')                 # Aileron area
c_a = Q_('0.3828 m')                # Aileron chord
delta_a = m.radians(30)     # [rad] max aileron deflection
delta_a *= Q_('rad')
delta_CL_max_a = 0.8267     # Max lift coeff difference due to aileron deflection
