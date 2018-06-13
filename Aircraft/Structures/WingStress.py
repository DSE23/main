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
from Geometry import Wing as GWing
import Wing
from Structures import Inertia
from Structures import Wing
from Aerodynamics import Wing as AWing
from Performance import Performance
import matplotlib.pyplot as plt



cl, cd, cm = AWing.computeloads()           #Load aerodynamic properties
n = 20                      #number of the devided sections
b = Geometry.Wing.b/2         #Wing span
b = b.magnitude * ureg.meter



z = Wing.z.magnitude * ureg.meter      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r.magnitude * ureg.meter      #root chord in m
rho = Performance.rho_c.magnitude * ureg("kg/(m**3)")         #cruise density
V = Performance.V_cruise.magnitude * ureg("m/s")        #cruise speed
zs = b - b/(n*2)     #subtract, zodat hij bij de eerste sectie op de helft begint
sectionlength = b/n
L_moment = Q_('0 kg * m ** 2 / s**2')
D_moment = Q_('0 kg * m ** 2 / s**2')
L = Q_('0 kg * m / s**2')
D = Q_('0 kg * m / s**2')
M = Q_('0 kg * m ** 2 / s**2')
dLlist = np.array([])
Llist = np.array([])
Dlist = np.array([])
zslist = np.array([])

while zs > z:                               #zs is measured is m from
    Areaofsection = sectionlength*Wing.length_chord(zs)
    dL = cl * 0.5 * rho * (V**2) * Areaofsection.magnitude * Q_("1 m**2")                #lift of the section
    dD = cd * 0.5 * rho * (V**2) * Areaofsection        #drag of the section
    dM = cm * 0.5 * rho * (V ** 2) * Areaofsection * Wing.length_chord(zs)      #moment of the section
    if zs < Geometry.Fuselage.b_f:
        dL = Q_('0 kg * m / s**2')
        dD = Q_('0 kg * m / s**2')
        dM = Q_('0 kg * m ** 2 / s**2')
    dL_moment = zs * dL                                 #moment produced by the lift on section
    dD_moment = zs * dD                                 #drag prduced by the drag on the section
    L = L + dL                      #Total lift for one wing
    D = D + dD                      #Total drag for one wing
    M = M + dM                      #Total moment for one wing
    L_moment = L_moment + dL_moment     #Total bending moment or
    D_moment = D_moment + dD_moment

    Llist = np.append(Llist, L)            #put the values in a list so we can plot them
    Dlist = np.append(Dlist, dD)
    zslist = np.append(zslist, abs(zs))

    zs = zs - sectionlength                 #Select other section for the next loop

Llist *= ureg("N/m")
Dlist *= ureg("N/m")
#print('L sum ', L)                  #print the values
#print('D sum ', D)
#print('M sum ', M)
#print('L_moment', L_moment)
#print('D_moment', D_moment)
#print("Llist=", Llist[0])
# plt.plot(zslist, Llist)
# plt.show()







#Material properties of the chosen material.
#Current chosen material:
#Carbon fiber reinforced carbon matrix composite (Vf:50%)
youngs_modulus = Q_("95 GPa")
yield_strength = Q_("18 MPa")
shear_modulus = Q_("36 GPa") #G





def Normal_stress_due_to_bending(cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx_wb*Inertia.Iyy_wb-Inertia.Ixy_wb**2
    inertia_term_1 = (Inertia.Iyy_wb*y-Inertia.Ixy_wb*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx_wb*cs-Inertia.Ixy_wb*y)/denominator_inertia_term
    sigma_zs = D_moment*inertia_term_1 + L_moment*inertia_term_2
    strain = sigma_zs /youngs_modulus 
    return strain #Gives the normal stress function for a given span zs, and x- and y- coordinate

# print('sigma_zs', Normal_stress_due_to_bending(0.15, Wing.airfoilordinate(Wing.c)))
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
def calc_moment_from_shear(qs, t_sk, t_sp, zs):
    Moment = 0
    step = 0.0001

    for x_c in np.arange(0, 1, 2*step):
        y_c_1 = Wing.airfoilordinate(x_c)
        y_c_2 = Wing.airfoilordinate(x_c+step)
        slope = (y_c_2 - y_c_1)/(step)
        Force_angle = -np.arctan(slope)
        perim = np.sqrt((step*Wing.length_chord(zs))**2 + ((y_c_2 - y_c_1)*Wing.length_chord(zs))**2)
        Loc_Force = perim * qs

        Fx = Loc_Force * np.cos(Force_angle)
        Fy = Loc_Force * np.sin(Force_angle)
        Moment += -Fx*(y_c_1 + y_c_2)*0.5*Wing.length_chord(zs)
        Moment += Fy*(x_c + 0.5*step)*Wing.length_chord(zs)

    return Moment

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
    qbase = Q_("0 N/m") #SHEAR IS NOT FINISHED
    return qs, qbase


def Torsion(zs, qbase):
    A_cell = Wing.Area_cell()
    length_skin = Wing.Area/Wing.ThSkin
    length_spar1 = Wing.airfoilordinate(Wing.ChSpar1)
    length_spar2 = Wing.airfoilordinate(Wing.ChSpar2)
    T = M #+ 2*(A_cell*qbase)
    const_tor = T/(4*A_cell**2*shear_modulus) #constant term in twist formula
    line_int_tor  = length_skin/Wing.ThSkin
    line_int_tor += length_spar1*Wing.length_chord(zs)/Wing.ThSpar1
    line_int_tor += length_spar2*Wing.length_chord(zs)/Wing.ThSpar2 #result from line integral from torsion formula
    twist_wb_tor_per_m  = const_tor*line_int_tor
    twist_wb_tor = twist_wb_tor_per_m*zs
    return twist_wb_tor
# print("twist =", Torsion(GWing.b/2,0))

# Wing deformation in X-direction
def deformation_x(zs):
    deformation_temp = Dlist[-1]/24*(zs)**4 #-Geometry.Fuselage.D_fus_max/2
    deformation_temp += -((Dlist[-1]-Dlist[0])/(GWing.b/2))/120*(zs)**5 #-Geometry.Fuselage.D_fus_max/2
    deformation_x = -deformation_temp/(youngs_modulus*Inertia.Ixx_wb)
    #deformation_x += D_moment/2*Geometry.Fuselage.D_fus_max**2/(youngs_modulus*Inertia.Ixx_wb)
    return deformation_x

#print("deformation_x=", deformation_x(GWing.b/2))


# Wing deformation in Y-direction
def deformation_y(zs):
    deformation_temp = Llist[-1]/24*(zs)**4 #-Geometry.Fuselage.D_fus_max/2
    deformation_temp += -((Llist[-1]-Llist[0])/(GWing.b/2))/120*(zs-Geometry.Fuselage.D_fus_max/2)**5 #-Geometry.Fuselage.D_fus_max/2
    deformation_y = -deformation_temp/(youngs_modulus*Inertia.Iyy_wb)
    #deformation_y += L_moment/2*Geometry.Fuselage.D_fus_max**2/(youngs_modulus*Inertia.Iyy_wb)
    return deformation_y.to(ureg("meter"))


#print("deformation_y=", deformation_y(GWing.b/2))



# with is like your try .. finally block in this case
with open('StrucVal.py', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

data[18] = 'youngs_modulus = Q_(\"' + str(youngs_modulus) + '\")\n'

# now change the 2nd line, note that you have to add a newline

# and write everything back
with open('StrucVal.py', 'w') as file:
    file.writelines(data)
