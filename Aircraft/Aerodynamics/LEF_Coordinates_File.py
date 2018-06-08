# -*- coding: utf-8 -*-
"""
Name: Leading Edge Flapped Airfoil Coordinate Calculator
Department: Aerodynamics
Last updated: 07/06/2018 09:47 by Emma
"""

# Imports

import math as m
import numpy as np
import scipy as sp
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

#import airfoil coordinates
airfoil_co = np.genfromtxt("../airfoil.dat")
#Variables to be defined by user

c_flap = 0.25# chord percentage covered by leading edge flap
c_hinge = 0.2# chord percentage of location hinge of leading edge flap
delta_flap = 15# [°] flap deflection in degrees
margin = 0.02 #tweak this to smoothen bottom part of airfoil
knuckle_smooth = 0.0735 #♠ tweak to smoothen knuckle, based on y-loc before flap


#Transforming data to SI units
delta_flap_rad = m.radians(delta_flap)

counter = 0
for coordinate in airfoil_co[:,0]:
    if coordinate <= c_flap:
        print('flap gets deflected')
        x = coordinate - c_hinge
        y = airfoil_co[counter, 1]
        #Euler transformation
        x_new =  c_hinge + (m.cos(delta_flap_rad) * x  -  m.sin(delta_flap_rad) * y)
        y_new = m.sin(delta_flap_rad) * x  +  m.cos(delta_flap_rad) * y
        if x_new < c_flap - margin :
            airfoil_co[counter,0] = x_new
            airfoil_co[counter,1] = y_new
        elif x_new >= c_flap - margin:
            print('------------Deleted stuff')
            airfoil_co = np.delete(airfoil_co, (counter), axis=0)
            counter = counter - 1
        if y_new > knuckle_smooth:
            print('smoothening knuckle')
            airfoil_co[counter,1] = knuckle_smooth
    else:
        print('Do nothing')
    print('counter', counter)
    counter = counter + 1
        
        
        
np.savetxt('new_airfoil.dat',airfoil_co)
