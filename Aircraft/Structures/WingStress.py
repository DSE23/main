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
from Structures import Wing
#from Aerodynamics import Wing as AWing
from Aerodynamics import Wing as AWing
from Performance import Performance


from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder


cl, cd, cm = AWing.computeloads()           #Load aerodynamic properties
n = 20                      #number of the devided sections
b = Geometry.Wing.b         #Wing span
z = Wing.z      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r      #root chord in m
rho = Performance.rho_c         #cruise density
V = Performance.V_cruise        #cruise speed
zs = b - b/(n*2)                   #subtract, zodat hij bij de eerste sectie op de helft begint
sectionlength = b/n
L_moment = Q_('0 kg * m ** 3 / s**2')
D_moment = Q_('0 kg * m ** 3 / s**2')


while zs > z:
    av_chord = (Wing.length_chord(z)+Wing.length_chord(b))/2        #average chord right from the crossection (m)
    spanleft = b - z
    Arealeft = spanleft*av_chord
    L = cl * 0.5 * rho * (V**2) * Arealeft
    D = cd * 0.5 * rho * (V**2) * Arealeft
    LoverSection = L * sectionlength                                             #Lift times the section length
    dL_moment = zs * LoverSection
    DoverSection = D * sectionlength
    dD_moment = zs * DoverSection
    L_moment = L_moment + dL_moment
    D_moment = D_moment + dD_moment

    zs = zs - sectionlength

print('L_moment', L_moment)
print('D_moment', D_moment)

M = cm * 0.5 * rho * (V**2) * Arealeft * Wing.length_chord(z)




#Material properties of the chosen material.
#Current chosen material:
#Carbon fiber reinforced carbon matrix composite (Vf:50%)
youngs_modulus = Q_("95 GPa")
yield_strength = Q_("18 MPa")
shear_modulus = Q_("36 GPa") #G





def Normal_stress_due_to_bending(zs, cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx_wb*Inertia.Iyy_wb-Inertia.Ixy_wb**2
    inertia_term_1 = (Inertia.Iyy_wb*y-Inertia.Ixy_wb*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx_wb*cs-Inertia.Ixy_wb*y)/denominator_inertia_term
    sigma_zs = D_moment*inertia_term_1 + L_moment*inertia_term_2
    return sigma_zs #Gives the normal stress function for a given span zs, and x- and y- coordinate

def Pure_torsion(zs):
    A_cell = Wing.Area_cell()
    length_skin = Wing.Area/Wing.ThSkin
    length_spar1 = Wing.airfoilordinate(Wing.ChSpar1)
    length_spar2 = Wing.airfoilordinate(Wing.ChSpar2)
    const_tor = M/(2*A_cell**2*shear_modulus) #constant term in twist formula
    line_int_tor = (length_skin*2/Wing.t_skin+length_spar1/Wing.ThSpar1+length_spar2/Wing.ThSpar2) #result from line integral from torsion formula
    twist_wb = const_tor*line_int_tor
    return twist_wb

