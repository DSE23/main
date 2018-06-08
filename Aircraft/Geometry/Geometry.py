# -*- coding: utf-8 -*-
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import math as m
import numpy as np

class Wing(object):
    
    S = Q_('11.74 m**2')                   # [m^2] Wing Surface
    A = 5.5                     # Aspect Ratio
    b = np.sqrt(S*A)             # [m] Wing Span
    taper = 0.45                # Taper ratio
    c_r = 2.015                 # Root chord
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
    
class Fuselage(object):
    
    l_f = Q_("6.22 m")                  # [m] Fuselage length
    D_fus_max = Q_("1.044 m")           # [m] Maximum fuselage diameter
    b_f = Q_("1.044 m")                 # [m] Fuselage width
    h_f = Q_("1.1 m")                   # [m] Fuselage height
    front_A = Q_("0.98 m**2")              # [m^2] Frontal area
    S_wet_f= Q_("14.94 m**2")              # [m^2] Fuselage wetted area
    
class H_tail(object):
    
    S = Q_("2.629 m**2")                 # [m^2] Horizontal tail surface
    A = 2.86                  # Aspect ratio
    b = np.sqrt(S*A)       # [m] Span horizontal tail
    taper = 0.529               # taper ratio
    c_r = Q_("1.329 m")                # [m] H-tail root chord
    c_t = c_r * taper         # [m] H-tail tip chord
    Sweep = Q_("0 deg")                 # [deg] Sweep H-tail
    MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))  # [m] Mean aerodynamic chord
    S_wet = S                 # [m^2] Wetted area
    S_e = Q_("1.3145 m**2")                # Elevator area
    c_e = Q_("0.50829 m")               # Elevator chord
    delta_e = m.radians(30)     # Max elevator deflection
    X_h = Q_("5.27 m")                  # [m] 0.25C location compared to the nose
    Z_h = Q_("0.55 m")                  # [m] Distance MAC_h and zero lift line wing
    
class V_tail(object):
    
    S = Q_("1.08 m**2")                                    # [m^2] V-tail surface area
    A = 0.99                                    # Aspect ratio
    b = np.sqrt(S*A)                             # [m] Span V-tail
    taper = 0.33                                # Taper ratio
    c_r = Q_("1.660 m")                                 # [m] V-tail root chord
    c_t = c_r * taper                           # [m] V-tail tip chord
    Sweep = Q_("0 deg")                                   # [deg] Sweep V-tail
    MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))  # [m] Mean aerodynamic chord
    S_wet = S                                   # [m^2] Wetted area
    S_r = Q_("0.5726 m**2")                                # Rudder area
    c_r = Q_("0.553 m")                                 # Rudder chord
    delta_r = m.radians(30)                     # Max rudder deflection
    X_v = Q_("5.70 m")                                  # [m] 0.25C location compared to the nose
    Z_v = Q_("0.52 m")                                  # [m] Distance MAC_h and zero lift line wing

    
class Landing_gear(object):
    
    lg_wheel_d = Q_("0.4445 m")                     # [m] Landing gear wheel diameter
    lg_wheel_w = Q_("0.16 m")                       # [m] Lg wheel width
