# -*- coding: utf-8 -*-
"""
Name: Slat Force
Department: Aerodynamics
Last updated: 13/06/2018 14:14 by Emma
"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
from aerocalc import std_atm

#Constants
slatlength = Q_(' m') #upper part of slat airfoil length
slatwidth = Q_('m') #width of the slats
Temp = Q_(' K') #Temperature at altitude you want to calc
density_0 = Q_('1.225 kg/m**3')
Temp0 = Q_('288.15 K')
g = Q_('9.80665 m / s**2')
R = Q_('287.05 m**2 * K / s**2')
lam = Q_('0.0065 K / m')
Velocity = Q_(' m/s') #aircraft velocity at which slats are deployed
p0 = Q_('101325 N/m**2')


# Slat force calculation
slatarea = slatlength * slatwidth
density = density_0 * (Temp / Temp0)**((g.magnitude/ lam.magnitude / R.magnitude) - 1)
p_stat = (Temp/Temp0)**(-g/(-.0065*R))*p0 #static pressure
P = (0.5 * density * Velocity**2 + p_stat) * slatarea #F = (dynamic + static pressure) * area
print('Slat force required =', P)