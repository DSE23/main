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
from Propulsion_and_system import Engine_mounts

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

y = 0
x = 0
z = 0                                         #the dreadful z
z *= Q_('meter')


#Boom areas section 1 (normal)
B_sec1 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2          #area for all booms


#Boom areas section 2 (cut out)
B_sec2 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2

#Boom areas section 3 (cone/taper)
def B_calc(x):
    b_f_taper = b_f80 / l_sec3 * x
    B_sec3 = t * (b_f_taper / 2) / 2 + t * (b_f_taper / 2) / 2
    return B_sec3, b_f_taper
B_sec3, b_f_taper = B_calc(x)

'''Inertia calculation'''

#Inertia section 1 (normal)
Izz_sec1, Iyy_sec1 = (b_f80/2)**2 * B_sec1 * 4

#Inertia section 2
Izz_sec2, Iyy_sec2 = (b_f80/2)**2 * B_sec2 * 4

#Inertia section 3
Izz_sec3, Iyy_sec3 = (b_f_taper/2)**2 * B_sec3 * 4

'''Loading force/moment calculation'''

#Loading section 1
My = x * Engine_mounts.f_y + Engine_mounts.m_y         #calculation of the resultant moment in y as variable of z
Mx = Engine_mounts.m_x         #calculation of the resultant moment in y as variable of z
Mz = x * Engine_mounts.f_z + Engine_mounts.f_z          #calculation of the resultant moment in y as variable of z

#Loading section 2




'''Bending stress calculations'''

#Bending section 1
sigma_x = My / Iyy_sec1 * z + Mz / Izz_sec1 * y         #bending stress

#Torsion section 2

q_x = Mx / (2 * b_f80)                                  #shear flow
shear_x = q_x / t

print(sigma_x)










