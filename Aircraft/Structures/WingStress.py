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
import Structures
from Structures import Inertia
from Structures import Wing

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder


cl = 1
cd = 0.05
cm = 0.05
b = Geometry.Wing.b         #Wing span
z = Structures.Wing.z      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r      #root chord in m

av_chord = (Wing.length_chord(z)+Wing.length_chord(b))/2        #average chord right from the crossection (m)
spanleft = b - z
print(av_chord)
print(spanleft)


#Material properties of the chosen material.
#Current chosen material:
#Carbon fiber reinforced carbon matrix composite (Vf:50%)
youngs_modulus = Q_("95 GPa")
yield_strength = Q_("18 MPa")
shear_modulus = Q_("36 GPa") #G



 

def Normal_stress_due_to_bending(zs, cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx*Inertia.Iyy-Inertia.Ixy**2
    inertia_term_1 = (Inertia.Iyy*y-Inertia.Ixy*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx*cs-Inertia.Ixy*y)/denominator_inertia_term
    sigma_zs = Mx*inertia_term_1 + My*inertia_term_2
    return sigma_zs #Gives the normal stress function for a given span zs, and x- and y- coordinate

def Pure_torsion(zs):
    A_cell = Wing.Area_cell()
    length_skin = Wing.Area/Wing.ThSkin
    length_spar1 = Wing.airfoilordinate(Wing.ChSpar1)
    length_spar2 = Wing.airfoilordinate(Wing.ChSpar2)
    const_tor = T/(2*A_cell**2*shear_modulus) #constant term in twist formula
    line_int_tor = (length_skin*2/Wing.t_skin+length_spar1/Wing.ThSpar1+length_spar2/Wing.ThSpar2) #result from line integral from torsion formula
    twist_wb = const_tor*line_int_tor
    return twist_wb