"""
Name: Geometry
Department: Geometry
Last updated: 08/06/2018 13:45 by Midas

In this file all of the geometric parameters of StefX should be present.
No other file should have these dimensions hard coded so all values are up to date

Example on how to use this file:
    from Geometry import Geometry

    S = Geometry.Wing.S
"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import math as m
import numpy as np
from Structures import StrucVal
from Propulsion_and_systems import Firewall
from Propulsion_and_systems import Engine
#from Propulsion_and_systems import Propeller

class Wing(object):
    
    S = Q_('11.74 m**2')                   # [m^2] Wing Surface
    A = 5.5                     # Aspect Ratio
    b = np.sqrt(S*A)             # [m] Wing Span
    taper = 0.45     # Taper ratio
    horn = Q_('0.0 m ')
    c_r = Q_("2.015 m")                 # Root chord
    c_t = c_r * taper           # Tip chord
    c_avg = (c_r + c_t)/2       #Average chord
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
    T_Cmax = 0.1513                             #Max thickness over chord
    S_wet = 2*S                   # Wetted wing area
    S_a = Q_('2.677 m**2')                 # Aileron area
    c_a = Q_('0.405 m')                # Aileron chord
    delta_a = Q_("30 deg")     # max aileron deflection
    delta_CL_max_a = 0.8267     # Max lift coeff difference due to aileron deflection
              
class Prop(object):
    Diameter = Q_('1.9 m')
    mass = Q_("30 kg")                   # Based on MT-propeller (4-bladed)

class Fuselage(object):
    
    l_f = Q_("6.22 m")                  # [m] Fuselage length
    D_fus_max = Q_("1.044 m")           # [m] Maximum fuselage diameter
    R_fus_firewall = Q_("0.50")         # [m] Fuselage radius at the Firewall
    b_f = Q_("1.044 m")                 # [m] Fuselage width
    h_f = Q_("1.1 m")                   # [m] Fuselage height
    front_A = Q_("0.98 m**2")              # [m^2] Frontal area
    S_wet_f= Q_("14.94 m**2")              # [m^2] Fuselage wetted area
    cabin_w = Q_("0.6 m")               # [m] cabin width (inside)
    A_max_canopy = Q_(" m**2")          # [m^2], coming from catia

class H_tail(object):
    
    S = Q_("1.8m**2")                 # [m^2] Horizontal tail surface
    A = 2.86                  # Aspect ratio
    b = np.sqrt(S*A)       # [m] Span horizontal tail
    taper = 0.529               # taper ratio
    c_r = 2 * (b/A)/(1+taper)                # [m] H-tail root chord
    c_t = c_r * taper         # [m] H-tail tip chord
    Sweep_25 = Q_("0 deg")                 # [deg] Sweep H-tail
    MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))  # [m] Mean aerodynamic chord
    S_wet = 2*S                 # [m^2] Wetted area
    S_e = Q_("1.3145 m**2")                # Elevator area
    ce_c = Q_("0.45")           # percentage elevator chord over local chord
    delta_e = Q_("30 deg")     # Max elevator deflection
    X_h = Q_("5.0 m")                  # [m] LE location compared to the nose
    Z_h = Q_("0.55 m")                  # [m] Distance MAC_h and zero lift line wing
    i_h = Q_("0 rad")                   #incidence angle ht
    Sweep_LE = (np.arctan(np.tan((Sweep_25))-(4/A) *
                       ((0-0.25)*(1-taper)/(1+taper))))
    
class V_tail(object):
    
    S = Q_("1.55 m**2")                                    # [m^2] V-tail surface area
    A = 0.99                                    # Aspect ratio
    b = np.sqrt(S*A)                             # [m] Span V-tail
    taper = 0.33                                # Taper ratio
    c_r = 2 * (b/A)/(1+taper)                                  # [m] V-tail root chord
    c_t = c_r * taper                           # [m] V-tail tip chord
    Sweep = Q_("0 deg")                                   # [deg] Sweep V-tail
    Sweep_LE = (np.arctan(np.tan((Sweep))-(4/A) *
                       ((0-0.25)*(1-taper)/(1+taper))))     #LE sweep angle
    MAC = c_r*(2/3)*((1+taper+taper**2)/(1+taper))  # [m] Mean aerodynamic chord
    S_wet = 2*S                                   # [m^2] Wetted area
    S_r = Q_("0.5726 m**2")                                # Rudder area
    cr_c = Q_("0.5")                            # percentage rudder chord over local chord
    delta_r = Q_("30 deg")                     # Max rudder deflection
    X_v = H_tail.X_h + Q_("0.43 m")                               # [m] 0.25C location compared to the nose
    Z_v = Q_("0.312 m")                                  # [m] Distance MAC_v and zero lift line wing
    t_c = 0.15                                  # [-] t/c V-tail
    
class Landing_gear(object):
    Prop_clear_req = Q_("0.23 m")                   # Required prop clearance CS23
    Prop_length = Prop.Diameter/2
    Z_mainlg = Prop_clear_req + Prop_length       # [m] Z_location bottom main L_G
    Tip_angle = Q_("12 deg")                        # Should be between 10 and 15
    X_mainlg = Firewall.xcg
    X_taillg = Fuselage.l_f
    Z_tailg = Z_mainlg - (X_mainlg-X_taillg) * np.tan(Tip_angle)
    Y_mainlg = Z_mainlg*(np.tan(Q_("35 deg")))
    lg_wheel_d = Q_("0.4445 m")                     # [m] Landing gear wheel diameter
    lg_wheel_w = Q_("0.16 m")                       # [m] Lg wheel width
    

class Masses(object):                    # !!!Structures should watch this!!!
    W_wing = StrucVal.Weightwing * 2     # Weight of the wing
    W_htail = Q_("20 kg")                # [kg] Mass of H_tail
    W_vtail = Q_("18 kg")                 # [kg] Mass of V_tail
    W_fus = Q_("82 kg")                  # [kg] Mass of Fuselage
    W_gear = Q_("60 kg")                 # [kg] Mass of landing gear
    W_engine = Engine.mass               # [kg] Mass of engine
    W_prop = Prop.mass                   # [kg] Mass of propellor
    W_fuelsys = Q_("10 kg")              # [kg] Mass of fuel system
    W_hydraulic = Q_("1 kg")             # [kg] Mass of hydraulics
    W_flightcontrol = Q_("20 kg")        # [kg] Mass of flight control
    W_avionics = Q_("17 kg")             # [kg] Mass of Avionics
    W_elecsys = Q_("46 kg")              # [kg] Mass of electronic systems
    W_lehld = Q_("16 kg")                # [kg] Mass of LE HLD's
    W_flaperons = Q_("14 kg")            # [kg] Mass of flaperons
    W_OEW = W_wing + W_htail + W_vtail + W_fus + W_gear + W_engine +\
            W_prop + W_fuelsys + W_hydraulic + W_flightcontrol +\
            W_avionics + W_elecsys + W_lehld + W_flaperons
    W_pilot = Q_("100 kg")              # [kg] Mass of pilot
    W_fuel = Q_("57 kg")                # [kg] Mass of fuel
    W_MTOW = W_OEW + W_pilot + W_fuel


class CG(object):
    # Locations in here are referencd to the nose in x and the crankshaft of the
    # engine in z, unless stated otherwise, For Z downward == positive, for x
    # is to the tail positive
    CG_wing_mac = 0.45                   # CG location of wing as percentage of MAC
    XLEMAC = Q_("1.1 m")               # LEMAC position
    CG_wing = CG_wing_mac*Wing.MAC + XLEMAC     # Wing CG position relative to nose
    X_wing = XLEMAC + 0.25 * Wing.MAC
    ZCG_wing = Q_("0.0 m")
    YCG_wing = Q_("0 m")
    CG_htail_mac = 0.5618
    CG_htail = H_tail.X_h + H_tail.MAC * CG_htail_mac   # H-tail cg relative to nose
    X_htail = H_tail.X_h + H_tail.MAC * 0.25            # 0.25c location of H_tail
    ZCG_htail = Q_("-0.094 m")
    YCG_vtail = Q_("0 m")
    CG_vtail_mac = 0.58799
    CG_vtail = V_tail.X_v + V_tail.MAC * CG_vtail_mac      # V-tail cg relative to nose
    ZCG_vtail = Q_("0.31m") - V_tail.b * 0.5
    YCG_vtail = Q_("0 m")
    CG_fus = Q_("2.16 m")                      # CG fuselage relative to nose !!!update!!!
    Z_fusorig = Q_("0.487 m")                  # Origin of fuselage in Z (lowest point)
    ZCG_fus = Q_("-0.0617m") #Inertia.ZCG_f                    # Complete fuselage Z-cg location
    YCG_fus = Q_("0 m")
    CG_lgear = 0.23 * Fuselage.l_f             # CG LG relative to nose !!!update!!!
    ZCG_lgear = Q_("0.69 m")
    CG_engine = Engine.xcg                     # CG of the engine relative to nose
    CG_prop = Q_("-20 cm")  # DUMMY VALUE, NOT KNOWN YET, negative because in front of datum
    ZCG_prop = Q_("0 m")
    ZCG_engine = Engine.zcg
    CG_fuelsys = Q_("1.18 m")                  # CG fuel system !!!update!!!
    ZCG_fuelsys = Q_("0 m")
    CG_hydraulics = Q_("1.115 m")              # CG hydraulics !!!update!!!
    ZCG_hydraulics = Q_("0 m")
    CG_flightcon = Q_("2.035 m")               # CG flight controls !!!update!!!
    ZCG_flightcon = Q_("0.21 m")
    CG_avionics = Q_("2.035 m")                # CG Avionics !!!update!!!
    ZCG_avionics = Q_("0 m")
    CG_elecsys = Q_("1.935 m")                 # CG electronic system !!!update!!!
    ZCG_elecsys = Q_("0 m")
    CG_lehld = XLEMAC                          # CG leading edge HLD's
    ZCG_lehld = Q_("0 m")
    CG_flaperons = XLEMAC + Wing.MAC           # CG Flaperons
    ZCG_flaperons = Q_("0 m")
    CG_pilot = Q_("2.235 m")                   # CG Pilot relative to nose
    ZCG_pilot = Q_("-0.22 m")
    CG_fuel = Q_("1.18 m")                     # CG fuel
    ZCG_fuel = Q_("0 m")
    CG_OEW = (Masses.W_wing * CG_wing + Masses.W_htail * CG_htail + Masses.W_vtail\
              * CG_vtail + Masses.W_fus * CG_fus + Masses.W_gear * CG_lgear\
              + Masses.W_engine * CG_engine + Masses.W_prop * CG_prop\
              + Masses.W_fuelsys * CG_fuelsys + Masses.W_hydraulic * CG_hydraulics\
              + Masses.W_elecsys * CG_elecsys + Masses.W_flightcontrol *\
              CG_flightcon + Masses.W_avionics * CG_avionics + Masses.W_lehld *\
              CG_lehld + Masses.W_flaperons * CG_flaperons)/Masses.W_OEW
    CG_wpilot = (CG_OEW * Masses.W_OEW + Masses.W_pilot * CG_pilot)/\
                (Masses.W_OEW+ Masses.W_pilot)
    CG_mtow = (CG_wpilot*(Masses.W_OEW+ Masses.W_pilot)+CG_fuel * Masses.W_fuel)\
              /(Masses.W_MTOW)
    ZCG_OEW = (Masses.W_wing * ZCG_wing + Masses.W_htail * ZCG_htail + Masses.W_vtail\
              * ZCG_vtail + Masses.W_fus * ZCG_fus + Masses.W_gear * ZCG_lgear\
              + Masses.W_engine * ZCG_engine + Masses.W_prop * ZCG_prop\
              + Masses.W_fuelsys * ZCG_fuelsys + Masses.W_hydraulic * ZCG_hydraulics\
              + Masses.W_elecsys * ZCG_elecsys + Masses.W_flightcontrol *\
              ZCG_flightcon + Masses.W_avionics * ZCG_avionics + Masses.W_lehld *\
              ZCG_lehld + Masses.W_flaperons * ZCG_flaperons)/Masses.W_OEW
    ZCG_wpilot = (ZCG_OEW * Masses.W_OEW + Masses.W_pilot * ZCG_pilot)/\
                (Masses.W_OEW+ Masses.W_pilot)
    ZCG_mtow = (ZCG_wpilot*(Masses.W_OEW+ Masses.W_pilot)+ZCG_fuel * Masses.W_fuel)\
              /(Masses.W_MTOW)
