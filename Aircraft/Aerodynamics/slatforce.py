# -*- coding: utf-8 -*-
"""
Name: Slat Force
Department: Aerodynamics
Last updated: 13/06/2018 14:14 by Emma
"""

#imports
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
import numpy as np
import math as m
from Geometry import Geometry

#Manual inputs
chordlength = .1  #determine part of chord used for slats
slatwidth = Geometry.Wing.b #width of the slats
h = Q_('100 m') #altitude of flight
Velocity = Q_(' 30 m/s') #aircraft velocity at which slats are deployed


#Calculate surface manin airfoil that will be slat (length in m)
airfoildata = np.genfromtxt("../airfoil.dat") #import airfoil
halfairfoil = airfoildata[0:81,:] # take upper half of airfoil coordinates
counter = -1
distance = np.array([])
for xco in halfairfoil[:,0]:
    counter = counter + 1
    if xco <= chordlength:
        print('xco', halfairfoil[counter,0])
        new_d = m.hypot(halfairfoil[counter,0] - halfairfoil[counter - 1,0],\
                        halfairfoil[counter,1] - halfairfoil[counter - 1,1])
        distance = np.append(distance,new_d)
        print(counter)
print('distance',distance)
slatlength = np.sum(distance) * Geometry.Wing.c_avg #upper part of slat airfoil length
print('slatlength',slatlength)


#ISA calculations
density_0 = Q_('1.225 kg/m**3')
Temp0 = Q_('288.15 K')
g = Q_('9.80665 m / s**2')
R = Q_('287.05 m**2 * K / s**2')
lam = Q_('0.0065 K / m')
p0 = Q_('101325 N/m**2')
Temp = Temp0 - lam * h
density = density_0 * (Temp / Temp0)**((g.magnitude/ lam.magnitude / R.magnitude) - 1)


# Slat force calculation
slatarea = slatlength * slatwidth
print('area',slatarea)
P = 0.5 * density * Velocity**2 * slatarea #F = pressure times area
print('Slat force required =', P)


#electric actuator size + stroke length
x_length = (chordlength - 0.008) * Geometry.Wing.c_avg.magnitude
print('x',x_length)
y_length = (Geometry.Wing.T_Cmax / 2 - 0.04) * Geometry.Wing.c_avg.magnitude
print('y',y_length)
angle = m.degrees(m.atan(y_length/x_length))
print(angle,' = angle')
stroke = m.sqrt(x_length**2 + y_length**2)
print('stroke length =', stroke)

x_attach = (0.184 - 0.2 * chordlength ) * Geometry.Wing.c_avg.magnitude
y_attach = 2/3 * Geometry.Wing.T_Cmax / 2 * Geometry.Wing.c_avg.magnitude
angle_actuator = m.degrees(m.atan(y_attach / x_attach))
print('actuator anagl', angle_actuator)
actuator = m.sqrt(x_attach**2 + y_attach**2)
print('actuator length',actuator)

