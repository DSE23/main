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


#Contracted slipstream diameter
V0 = 30 #free stream velocity
data = Propeller.Thrustcalc(V0)
Dia = Geometry.Prop.Diameter.magnitude
Sref = Geometry.Wing.S.magnitude
CT = data[1].magnitude / 0.5 / 1.225 / V0**2 / Sref

DeltaV = V0 * (m.sqrt(1 + CT * (Sref / (Dia**2 / 4 * m.pi)) - 1))

# D_con = Dia * m.sqrt((V0 + DeltaV / 2)/(V0 + DeltaV))