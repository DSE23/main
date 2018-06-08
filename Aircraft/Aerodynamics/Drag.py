# -*- coding: utf-8 -*-
"""
Name: DATCOM Drag
Department: Aerodynamics
Last updated: 08/06/2018 16:53 by Emma
"""

#Imports

import math as m
import numpy as np
from Geometry import Geometry


#Input variables, to be filled in by user


Temp = # Temperature in Kelvin
Temp *= Q_('K')


#Standard values Inputs


density_0 = Q_('1.225 kg/m**3')
Temp0 = Q_('288.15 K')
g = Q_('9.80665 m / s**2')
R = Q_('287.05 m**2 * K / s**2')
lam = Q_('0.0065 K / m')


#Input variables, linked to aircraft geometry files etc.


chord = Geometry.Wing.c_avg #average chord



#Reynolds number calculation


viscosity = 1.458e**-6 * Temp**1.5 * (1 / (Temp + 110.4))
viscosity *= Q_('N * s / m**2')
density = density_0 * (Temp / Temp0)**((g/ lam / R) - 1)
Reynolds = density * velocity * chord / viscosity