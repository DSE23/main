# -*- coding: utf-8 -*-
"""
Name: Joukowsky Airfoil
Department: Aerodynamics
Last updated: 22/06/2018 15:03 by Emma
"""
import numpy as np
import math as m
thickness_range = np.linspace(0.1,0.15,20)
C_range = np.linspace(0,1,20)
chord = 1 #dummy variable!
for eps in thickness_range:
    C = 4 * chord / (3 + 2 * eps + (1/(1 + 2 * eps)))
    radius = C / 4 *( 1 + eps)
    zmax = 0
    for theta in np.linspace(m.pi/2,m.pi,30):
        z = radius * m.sin(theta) * (1 - (C**2 / 16)/((radius * m.cos(theta) -\
                                          eps * C /4)**2 + radius ** 2 *\
                                          (m.sin(theta)**2)))
        if z > zmax:
            zmax = z
    if zmax <= 0.076:
        print('eps', eps, 'zmax', zmax)
    else:
        break