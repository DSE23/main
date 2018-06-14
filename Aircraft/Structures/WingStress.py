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
n = 40                      #number of the devided sections
b = Geometry.Wing.b/2         #Wing span
b = b.magnitude * ureg.meter



z = Wing.z.magnitude * ureg.meter      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r.magnitude * ureg.meter      #root chord in m
rho = Performance.rho_c.magnitude * ureg("kg/(m**3)")         #cruise density
V = Performance.V_cruise.magnitude * ureg("m/s")        #cruise speed

dL_moment = Q_('0 kg * m ** 2 / s**2')
dD_moment = Q_('0 kg * m**2/s**2')
L_moment = Q_('0 kg * m ** 2 / s**2')
D_moment = Q_('0 kg * m ** 2 / s**2')

L = Q_('0 kg * m / s**2')
D = Q_('0 kg * m / s**2')
M = Q_('0 kg * m ** 2 / s**2')

Lmomentlist = np.array([])
Dmomentlist = np.array([])
dLlist = np.array([])
dDlist = np.array([])
Llist = np.array([])
Dlist = np.array([])

Sectioncenters = np.array([])

sectionlength = b/n

zs = b - (sectionlength/2)

while zs > z:                            #zs is measured is m from

    Areaofsection = sectionlength*Wing.length_chord(zs)

    Sectioncenters = np.append(Sectioncenters, zs)

    '''Lift drag and moment for the section'''
    dL = cl * 0.5 * rho * (V**2) * Areaofsection.magnitude * Q_("m**2")                #lift of the section
    dD = cd * 0.5 * rho * (V**2) * Areaofsection        #drag of the section
    dM = cm * 0.5 * rho * (V ** 2) * Areaofsection * Wing.length_chord(zs)      #moment of the section

    if zs < Geometry.Fuselage.b_f:
        dL = Q_('0 kg * m / s**2')
        dD = Q_('0 kg * m / s**2')
        dM = Q_('0 kg * m ** 2 / s**2')

    dLlist = np.append(dLlist, dL)
    dDlist = np.append(dDlist, dD)

    '''Total lift, drag and moment for the wing'''
    L = L + dL  # Total lift for one wing
    D = D + dD  # Total drag for one wing
    M = M + dM  # Total moment for one wing

    Llist = np.append(Llist, L)  # put the values in a list so we can plot them
    Dlist = np.append(Dlist, D)

    zs = zs - sectionlength  # Select other section for the next loop

for i in np.arange(0, len(dLlist)):
    arm = Sectioncenters[i] - z.magnitude
    dLmoment = arm * dLlist[i]
    dDmoment = arm * dDlist[i]
    L_moment = L_moment.magnitude + dLmoment
    D_moment = D_moment.magnitude + dDmoment
    Lmomentlist = np.append(Lmomentlist, L_moment)
    Dmomentlist = np.append(Dmomentlist, D_moment)
    L_moment *= Q_('kg * m**2/s**2')
    D_moment *= Q_('kg * m**2/s**2')


'''For the 20G manoeuver'''
MTOW = Geometry.Masses.W_MTOW
Max_20G_N = MTOW * 9.81 * 20
Tot_L = 2 * L
if Tot_L.magnitude > 0.:
    fac_20G = Max_20G_N / Tot_L
    fac_20G = fac_20G.magnitude
else:
    fac_20G = 0

L_moment = L_moment * fac_20G
D_moment = D_moment * fac_20G
L = L * fac_20G
D = D * fac_20G
M = M * fac_20G


Llist *= ureg("N/m")
Dlist *= ureg("N/m")
#print('L sum ', L)                  #print the values
#print('D sum ', D)
#print('M sum ', M)
#print('L_moment', L_moment)
#print('D_moment', D_moment)
#print("Llist=", Llist[0])
# plt.plot(Sectioncenters, Lmomentlist)
# plt.show()







#Material properties of the chosen material.
#Current chosen material:
#Carbon fiber reinforced carbon matrix composite (Vf:50%)
youngs_modulus = Q_("95 GPa") #E
yield_strength = Q_("23 MPa") #tensile
compr_strength = Q_("247 MPa") #compression
shear_modulus = Q_("36 GPa") #G





def Normal_stress_due_to_bending(cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx_wb*Inertia.Iyy_wb-Inertia.Ixy_wb**2
    inertia_term_1 = (Inertia.Iyy_wb*y-Inertia.Ixy_wb*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx_wb*cs-Inertia.Ixy_wb*y)/denominator_inertia_term
    sigma_zs = D_moment*inertia_term_1 + L_moment*inertia_term_2
    strain = sigma_zs /youngs_modulus
    return sigma_zs #Gives the normal stress function for a given span zs, and x- and y- coordinate

# print('sigma_zs', Normal_stress_due_to_bending(0.15, Wing.airfoilordinate(Wing.c)))
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED
#SHEAR IS NOT FINISHED


def Shear_wb(zs, L):
    #section 01
    n = 100
    ds = Wing.HSpar1/n
    qs1 = np.array([])
    s1 = np.array([])
    qs1 = np.append(qs1, 0)
    s1 = np.append(s1, 0)
    s = 0
    for i in range(n):
        s = s + ds
        s1 = np.append(s1, s)
        qs = s**2*Wing.ThSpar1*(-L)/Inertia.Ixx_wb
        qs1 = np.append(qs1, qs)
    section01at1 = np.sum(qs1)
    #section12
    n = 100 #number of sections
    ds = (Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)/n)
    qs2 = np.array([])
    s2 = np.array([])
    qs2 = np.append(qs2, section01at1)
    s2 = np.append(s2, 0)
    line_int_skin_wb = section01at1
    s = 0
    for i in range(n):
        s = s + ds
        s2 = np.append(s2, s)
        qs = s * Wing.get_xy_from_perim(s/Wing.length_chord(zs))[1]*Wing.length_chord(zs)*Wing.ThSkin*(-L)/Inertia.Ixx_wb
        qs2 = np.append(qs2, qs)
    section12at2 =  np.sum(qs2)
    #section23
    n = 100
    ds = Wing.HSpar2/n
    qs3 = np.array([])
    s3 = np.array([])
    qs3 = np.append(qs3, section12at2)
    s3 = np.append(s3, 0)
    s = 0
    for i in range(n):
        s = s + ds
        s3 = np.append(s3, s)
        qs = -s**2*Wing.ThSpar2*(-L)/Inertia.Ixx_wb
        qs3 = np.append(qs3, qs)
    qbase = 2 #*(Moment """+ hier moet nog iets bij""")/Wing.Area_cell()
    qs1 *= ureg("N/m")
    qs2 *= ureg("N/m")
    qs3 *= ureg("N/m")
    s1 *= ureg("m")
    s2 *= ureg("m")
    s3 *= ureg("m")
    shear_arrays = np.array([[qs1, s1], [qs2, s2], [qs3, s3]])
    return qs, qbase, shear_arrays

qs, qbase, shear_arrays = Shear_wb(Wing.z, L)
print("q1=", shear_arrays[0,0])
print("s1=", shear_arrays[0,1])
print("q2=", shear_arrays[1,0])
print("s2=", shear_arrays[1,1])
print("q3=", shear_arrays[2,0])
print("s3=", shear_arrays[2,1])

def calc_moment_from_shear(zs, shear_arrays=Shear_wb(Wing.z, L)[2]):
    q_1 = shear_arrays[0, 0]
    s_1 = shear_arrays[0, 1]
    q_2 = shear_arrays[1, 0]
    s_2 = shear_arrays[1, 1]
    q_3 = shear_arrays[2, 0]
    s_3 = shear_arrays[2, 1]
    Moment = 0
    print(zs)
    if (len(q_2) != len(s_2)):
        raise ValueError("ERROR, ARRAY LENGTH NOT EQUAL!")
    for i in range(0,len(s_2)-1,2):
        ds = s_2[i+1] - s_2[i]
        q_loc = (q_2[i] + q_2[i+1])*0.5
        #print(s_2[i])
        x_loc_1, y_loc_1 = Wing.get_xy_from_perim(s_2[i]/Wing.length_chord(zs), Wing.ChSpar1)
        x_loc_2, y_loc_2 = Wing.get_xy_from_perim(s_2[i+1]/Wing.length_chord(zs), Wing.ChSpar1)
        x_loc_1 *= Wing.length_chord(zs)
        y_loc_1 *= Wing.length_chord(zs)
        x_loc_2 *= Wing.length_chord(zs)
        y_loc_2 *= Wing.length_chord(zs)
        slope = (y_loc_2 - y_loc_1)/(x_loc_2 - x_loc_1)
        force_angle = -np.arctan(slope)
        print("q=", q_loc)
        print("ds=", ds)
        print("Angle=", force_angle)
        Fx = q_loc*ds*np.cos(force_angle)
        Fy = q_loc*ds*np.sin(force_angle)
        Moment += -Fx*(y_loc_1+y_loc_2)*0.5
        Moment += Fy*(x_loc_1+x_loc_2)*0.5

    return Moment

print(calc_moment_from_shear(z, shear_arrays))

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


#Tsia-Wu Failure criterion
def Tsia_Wu(sigma_zs, shearforce):
    F11=1/(yield_strength*compr_strength)
    F22 = F11
    F12 = -1/2*np.sqrt(F11*F22)
    F1 = 1/(yield_strength)-1/(compr_strength)
    F2 = 1/(yield_strength)-1/(compr_strength)
    #F44 = 1
    #F66 = 1
    F = F11 *F22*F12*F1*F2 #klopt niet
    return sigma_zs




# with is like your try .. finally block in this case
with open('StrucVal.py', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

data[18] = 'youngs_modulus = Q_(\"' + str(youngs_modulus) + '\")\n'

# now change the 2nd line, note that you have to add a newline

# and write everything back
with open('StrucVal.py', 'w') as file:
    file.writelines(data)

