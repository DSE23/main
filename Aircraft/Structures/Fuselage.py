"""                  
Name: Fuselage
Department: Structures
Last updated: 05/06/2018 12:45 by Midas
"""

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from scipy import optimize
from scipy import interpolate
from Geometry import Geometry

from Misc import ureg, Q_ # Imports the unit registry from the Misc folder

b_f = Geometry.Fuselage.b_f                   #diameter of fuselage at the fire wall
b_f80 = b_f*0.8                               #design radius
t = Q_('0.002 m')                             #thickness of fuselage skin
moter_w = Q_('180 kg')                          #Resultant weight of the motor
g = Q_('9.81 m / s**2')                         #The gravity acceleration
l_fus = Geometry.Fuselage.l_f                   #length of the fuselage in m
l_sec1 = Q_('1.0 m')                          #length of section 1 (normal)
l_sec2 = Q_('1.0 m')                          #length of section 2 (cut out)
l_sec3 = Q_('2.0 m')                          #length of section 3 (taper)

x =
z = 0                                         #the dreadful z
z *= Q_('meter')


#Boom areas section 1 (normal)
B_sec1 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2          #area for all booms


#Boom areas section 2 (cut out)
B_sec2 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2

#Boom areas section 3 (cone/taper)
def B_calc(z):
    b_f_taper = b_f80 / l_sec3 * z
    B_sec3 = t * (b_f_taper / 2) / 2 + t * (b_f_taper / 2) / 2
    return B_sec3, b_f_taper
B_sec3, b_f_taper = B_calc(z)

'''Inertia calculation'''

#Inertia section 1 (normal)
Ixx_sec1, Iyy_sec1 = (b_f80/2)**2 * B_sec1 * 4

#Inertia section 2
Ixx_sec2, Iyy_sec2 = (b_f80/2)**2 * B_sec2 * 4

#Inertia section 3
Ixx_sec3, Iyy_sec3 = (b_f_taper/2)**2 * B_sec3 * 4

'''Loading force/moment calculation'''

#Loading section 1
My = (l_sec1 - z) * (g * moter_w)           #calculation of the resultant moment in y as variable of z
Mx = 0



'''Bending stress calculations'''

#Bending section 1
sigma_z = My / Iyy_sec1 * x + Mx / Ixx_sec1 * y





