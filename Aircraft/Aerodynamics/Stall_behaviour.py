# -*- coding: utf-8 -*-
"""
Name: Stall Behaviour
Department: Aerodynamics
Last updated: 18/06/2018 15:03 by Emma
"""

#imports
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
import math as m
import numpy as np
import scipy as sp
from Geometry import Geometry
from pint import UnitRegistry
import Aeroprops
from Propulsion_and_systems import Propeller


#import propeller data
axial_v = np.genfromtxt('../DataWalong.txt')
wlist = np.array([])
range_for = len(axial_v[1,:])
for line in range(0,range_for):
    w = np.average(axial_v[line,:])
    wlist = np.append(wlist, w)
    

#Aircraft Geometry imports
width_fus = Geometry.Fuselage.D_fus_max.magnitude / 2
LE_wing = 1
root_c = Geometry.Wing.c_r.magnitude
aileron_c = Geometry.Wing.c_a.magnitude
LE_ht = Geometry.H_tail.X_h.magnitude
root_c_ht = Geometry.H_tail.c_r.magnitude
elevator_c_ratio = Geometry.H_tail.ce_c.magnitude
taper_ht = Geometry.H_tail.taper
span_ht = Geometry.H_tail.b.magnitude
LE_vt = Geometry.V_tail.X_v.magnitude
root_c_vt = Geometry.V_tail.c_r.magnitude
taper_vt = Geometry.V_tail.taper
rudder_c_ratio = Geometry.V_tail.cr_c.magnitude
span_vt = Geometry.V_tail.b.magnitude
#Contracted slipstream diameter
V0 = 30                                                                 #free stream velocity
Dia = Geometry.Prop.Diameter.magnitude                                  #diameter of the prop
Radius_p = Dia / 2
#ask Gijs which a to select for which velocity
a = wlist[1] / V0
l = Geometry.Fuselage.l_f
for x in np.linspace(0,l,200):
    Radius_tube = Radius_p * m.sqrt((1 + a) / (1 + a * (1 + (x/(m.sqrt(Radius_p ** 2 + x **2))))))
    velocity = a*V0 + V0
    #print(Radius_tube)
    
    if LE_wing <= x <= LE_wing + 0.03:
        area_ail = (Radius_tube - width_fus) * aileron_c * 2
        velocity_ail = velocity
        print('aileron area in flow:', area_ail, 'velocity of flow :', velocity_ail)
        
    elif LE_ht <= x <= LE_ht + 0.03:
        lam_local = taper_ht + (1 - taper_ht)/span_ht * 2 * (span_ht/2 - Radius_tube)
        chord_l = lam_local * root_c_ht * elevator_c_ratio
        area_elevator = 2 * (chord_l + root_c_ht * elevator_c_ratio) * Radius_tube / 2
        velocity_ele = velocity
        print('elevator area in flow:', area_elevator, 'velocity of flow :', velocity_ele)
        
    elif LE_vt <= x <= LE_vt + 0.03:
        lam_local = taper_vt + (1 - taper_vt)/span_vt * (span_vt - Radius_tube)
        chord_l = lam_local * root_c_vt * rudder_c_ratio
        area_rudder = (chord_l + root_c_vt * rudder_c_ratio) * Radius_tube / 2
        if area_rudder > Geometry.V_tail.S_r.magnitude:
            area_rudder = Geometry.V_tail.S_r.magnitude
        velocity_rud = velocity
        print('rudder area in flow:', area_rudder, 'velocity of flow :', velocity_rud)




