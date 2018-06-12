"""                  
Name: Wing_Stress_Calculations 
Department: Structures 
Last updated: 08/06/2018 12:38 by Boris 
"""

## Forces
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
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
import matplotlib.pyplot as plt



cl, cd, cm = AWing.computeloads()           #Load aerodynamic properties
n = 20                      #number of the devided sections
b = Geometry.Wing.b/2         #Wing span
z = Wing.z      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r      #root chord in m
rho = Performance.rho_c         #cruise density
V = Performance.V_cruise        #cruise speed
zs = b - b/(n*2)     #subtract, zodat hij bij de eerste sectie op de helft begint
sectionlength = b/n
L_moment = Q_('0 kg * m ** 2 / s**2')
D_moment = Q_('0 kg * m ** 2 / s**2')
L = Q_('0 kg * m / s**2')
D = Q_('0 kg * m / s**2')
M = Q_('0 kg * m ** 2 / s**2')
Llist = np.array([])
zslist = np.array([])

while zs > z:                               #zs is measured is m from
    Areaofsection = sectionlength*Wing.length_chord(zs)
    dL = cl * 0.5 * rho * (V**2) * Areaofsection                #lift of the section
    dD = cd * 0.5 * rho * (V**2) * Areaofsection        #drag of the section
    dM = cm * 0.5 * rho * (V ** 2) * Areaofsection * Wing.length_chord(zs)      #moment of the section
    dL_moment = zs * dL                                 #moment produced by the lift on section
    dD_moment = zs * dD                                 #drag prduced by the drag on the section
    L = L + dL                      #Total lift for one wing
    D = D + dD                      #Total drag for one wing
    M = M + dM                      #Total moment for one wing
    L_moment = L_moment + dL_moment     #Total bending moment or
    D_moment = D_moment + dD_moment

    Llist = np.append(Llist, dL)            #put the values in a list so we can plot them
    zslist = np.append(zslist, abs(zs))

    zs = zs - sectionlength                 #Select other section for the next loop

print('L sum ', L)                  #print the values
print('D sum ', D)
print('M sum ', M)
print('L_moment', L_moment)
print('D_moment', D_moment)

# plt.plot(zslist, Llist)
# plt.show()







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


def Shear_wb(zs):
    #section 01
    section01at1 = Wing.ThSpar1*Wing.HSpar1**2
    #section12
    n = 100 #number of sections
    ds = (Wing.arclength/n)
    s = 0
    line_int_skin_wb = section01at1
    for i in range(n):
        s = s + ds
        dline_int_skin_wb = s * Inertia.get_y_for_perimeter(x)
        line_int_skin_wb += dline_int_skin_wb
    section12at2 = Wing.ThSkin * line_int_skin_wb
    #section23
    section23at3 = section12at2 - Wing.ThSpar2*Wing.HSpar2**2
    qs = -L/Inertia.Ixx_wb*(section01at1+section12at2+section23at3)
    qbase = 1
    return qs, qbase


def Torsion(qbase):
    A_cell = Wing.Area_cell()
    length_skin = Wing.Area/Wing.ThSkin
    length_spar1 = Wing.airfoilordinate(Wing.ChSpar1)
    length_spar2 = Wing.airfoilordinate(Wing.ChSpar2)
    T = M + 2*A_cell*qbase
    const_tor = T/(2*A_cell**2*shear_modulus) #constant term in twist formula
    line_int_tor = (length_skin*2/Wing.t_skin+length_spar1/Wing.ThSpar1+length_spar2/Wing.ThSpar2) #result from line integral from torsion formula
    twist_wb_tor = const_tor*line_int_tor
    return twist_wb_tor

# Wing deformation in X-direction
def deformatio_x(zs):
    deformation_temp = drag_at_root/24*(zs-widthfuselage)^4
    deformation_temp += drag_slope/120*(zs-widthfuselage)^5
    deformation_x = 1/(youngs_modulus*Inertia.Ixx_wb)*deformation_temp
    deformation_x += L_moment/2*widthfuselage
    return deformation_x


def deformatio_y(zs):
    deformation_temp = lift_at_root/24*(zs-Geometry.D_fus_max/2)^4
    deformation_temp += lift_slope/120*(zs-Geometry.D_fus_max/2)^5
    deformation_y = 1/(youngs_modulus*Inertia.Ixx_wb)*deformation_temp
    deformation_y += D_moment/2*Geometry.D_fus_max/2
    return deformation_y
