"""                  
Name: Wing_Stress_Calculations 
Department: Structures 
Last updated: 08/06/2018 12:38 by Boris 
"""

## Forces
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

import numpy as np
from scipy import interpolate
import math as m
from Geometry import Geometry
import Wing
from Structures import Inertia
#from Aerodynamics import Wing as AWing
from Performance import Performance

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder


cl_max = 1 #AWing.CL_max

cd = 0.05
cm = 0.05
b = Geometry.Wing.b         #Wing span
z = Wing.z      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r      #root chord in m
rho = Performance.rho_c
V = Performance.V_cruise

av_chord = (Wing.length_chord(z)+Wing.length_chord(b))/2        #average chord right from the crossection (m)
spanleft = b - z
Arealeft = spanleft*av_chord
L = cl_max * 0.5 * rho * (V**2) * Arealeft

print(L)


def Normal_stress_due_to_bending(zs, cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx*Inertia.Iyy-Inertia.Ixy**2
    inertia_term_1 = (Inertia.Iyy*y-Inertia.Ixy*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx*cs-Inertia.Ixy*y)/denominator_inertia_term
    sigma_zs = Mx*inertia_term_1 + My*inertia_term_2
    return sigma_zs #Gives the normal stress function for a given span zs, and x- and y- coordinate


