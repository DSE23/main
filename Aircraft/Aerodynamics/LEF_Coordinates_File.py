# -*- coding: utf-8 -*-
"""
Name: Leading Edge Flapped Airfoil Coordinate Calculator
Department: Aerodynamics
Last updated: 06/06/2018 12:01 by Emma
"""

# Imports

import math as m
import numpy as np
import scipy as sp
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

#import airfoil coordinates
airfoil_co = np.genfromtext('../airfoil.dat')
#Variables to be defined by user

c_flap = # chord percentage covered by leading edge flap
c_hinge = # chord percentage of location hinge of leading edge flap
delta_flap = # [°] flap deflection in degrees


#Transforming data to SI units
delta_flap_rad = m.radians(delta_flap)

counter = 0
for coordinate in airfoil_co[:,0]:
    if coordinate <= c_flap:
        x = coordinate - c_hinge
        y = airfoil_co[counter, 1]
        #Euler transformation
        x_new = m.cos(delta_flap_rad) * x  -  m.sin(delta_flap_rad) * y
        y_new = m.sin(delta_flap_rad) * x  +  m.cos(delta_flap_rad) * y
        if x_new <= c_flap:
            airfoil_co[counter,0] = x_new
            airfoil_co[counter,1] = y_new
        elif x_new > c_flap:
            np.delete(airfoil_co, counter, 0)
            np.delete(airfoil_co, counter, 1)
    counter = counter + 1
        
        
        
np.savetxt('new_airfoil.dat')
