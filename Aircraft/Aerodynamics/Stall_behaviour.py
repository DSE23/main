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

axial_v = np.genfromtxt('../DataWalong')
wlist = np.array([])
for line in len(axial_v[1,:]):
    w = np.average(axial_v[line,:])
    wlist = np.append(alist, w)
#Contracted slipstream diameter
V0 = 30                                                                 #free stream velocity
Dia = Geometry.Prop.Diameter.magnitude                                  #diameter of the prop
Radius_p = Dia / 2
#ask Gijs which a to select for which velocity
a = w / V0
l = Geometry.Fuselage.l_f
for x in np.linspace(0,l,200):
    Radius_tube = Radius_p * m.sqrt((1 + a) / (1 + a * (1 + (x/(m.sqrt(Radius_p ** 2 + x **2))))))
    velocity = a*V0 + V0
