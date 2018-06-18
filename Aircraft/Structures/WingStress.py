"""                  
Name: Wing_Stress_Calculations 
Department: Structures 
Last updated: 15/06/2018 11:17 by Midas
"""

## Forces
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import numpy as np
from scipy import interpolate
import math as m
from Geometry import Geometry
# from Geometry import Wing as GWing
# import Wing
from Structures import Inertia
from Structures import Wing
from Aerodynamics import Wing as AWing
from Performance import Performance
import matplotlib.pyplot as plt
import time


cl, cd, cm = AWing.computeloads()           #Load aerodynamic properties
n = 10                      #number of the devided sections
b = Wing.s         #Wing span
b = b.magnitude * ureg.meter
z = Wing.z

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

sectionlength = (b.magnitude-0*Geometry.Fuselage.b_f.magnitude)/n*ureg.meter

zs = b - (sectionlength/2)

while zs.magnitude > 0:  # zs is measured is m from
    Areaofsection = sectionlength * Wing.length_chord(zs)

    Sectioncenters = np.append(Sectioncenters, zs)

    '''Lift drag and moment for the section'''
    dL = cl * 0.5 * rho * (V ** 2) * Areaofsection.magnitude * Q_("m**2")  # lift of the section

    if zs < Geometry.Fuselage.b_f * 0:
        dL = Q_('0 kg * m / s**2')

    '''Total lift, drag and moment for the wing'''
    L = L + dL  # Total lift for one wing

    zs = zs - sectionlength  # Select other section for the next loop

totallift = L


def computeloads(z):
    Sectioncenters = np.array([])
    dLlist = np.array([])
    dDlist = np.array([])
    L = 0
    D = 0
    M = 0
    L_moment = 0
    D_moment = 0
    Llist = np.array([])
    Dlist = np.array([])
    Lmomentlist = np.array([])
    Dmomentlist = np.array([])
    L_momentlist = np.array([])
    zs = b - (sectionlength / 2)
    while zs > z:                            #zs is measured is m from

        Areaofsection = sectionlength*Wing.length_chord(zs)

        Sectioncenters = np.append(Sectioncenters, zs)

        '''Lift drag and moment for the section'''
        dL = (cl * 0.5 * rho * (V**2) * Areaofsection).magnitude              #lift of the section
        dD = (cd * 0.5 * rho * (V**2) * Areaofsection).magnitude     #drag of the section
        dM = (cm * 0.5 * rho * (V ** 2) * Areaofsection * Wing.length_chord(zs)).magnitude     #moment of the section
        dL *= Q_('kg * m / s**2')
        dD *= Q_('kg * m / s**2')
        dM *= Q_('kg * m**2 / s**2')
        
        if zs < Geometry.Fuselage.b_f*0:
            dL = 0
            dD = 0
            dM = 0

        dLlist = np.append(dLlist, dL)
        dDlist = np.append(dDlist, dD)

        '''Total lift, drag and moment for the wing'''
        L = L + dL  # Total lift for one wing
        D = D + dD  # Total drag for one wing
        M = M + dM  # Total moment for one wing

        Llist = np.append(Llist, L)  # put the values in a list so we can plot them
        Dlist = np.append(Dlist, D)

        zs = zs - sectionlength  # Select other section for the next loop

    for i in range(0, len(Sectioncenters)):
        arm = (Sectioncenters[i] - z.magnitude)
        dLmoment = (arm * dLlist[i])
        dDmoment = (arm * dDlist[i])
        L_moment = L_moment + dLmoment
        D_moment = D_moment + dDmoment
        Lmomentlist = np.append(Lmomentlist, L_moment)
        Dmomentlist = np.append(Dmomentlist, D_moment)

    '''For the 20G manoeuver'''
    MTOW = Geometry.Masses.W_MTOW
    Max_20G_N = MTOW * 9.81 * 20
    Tot_L = 2 * totallift
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

    L *= Q_('kg * m / s**2')
    D *= Q_('kg * m / s**2')
    M *= Q_('kg * m ** 2 / s**2')
    L_moment *= Q_('kg * m ** 2 / s**2')
    D_moment *= Q_('kg * m ** 2 / s**2')

    return L, D, M, L_moment, D_moment, dL, dD, dM

L, D, M, L_moment, D_moment, dL, dD, dM = computeloads(z)
print('L', L)
print('D', D)
print('M', M)
#
#
# Llist *= ureg("N/m")
# Dlist *= ureg("N/m")
# # print('L sum ', L)                  #print the values
# # print('D sum ', D)
# # print('M sum ', M)
# # print('D_moment', D_moment)
# # print("Llist=", Llist[0])
# plt.plot(zlist, L_momentlist)
# plt.show()







#Material properties of the chosen material.
#Current chosen material:
#Carbon fiber reinforced carbon matrix composite (Vf:50%)
youngs_modulus = Q_("95 GPa")  #E
yield_strength = Q_("23 MPa")  #tensile
compr_strength = Q_("247 MPa") #compression
shear_modulus = Q_("36 GPa")   #G
poisson = 0.31                 # maximum 0.33
tau_max = Q_("35 MPa")





def Normal_stress_due_to_bending(cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx_wb*Inertia.Iyy_wb-Inertia.Ixy_wb**2
    inertia_term_1 = (Inertia.Iyy_wb*y-Inertia.Ixy_wb*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx_wb*cs-Inertia.Ixy_wb*y)/denominator_inertia_term
    sigma_zs = D_moment*inertia_term_1 + L_moment*inertia_term_2
    strain = sigma_zs /youngs_modulus
    return sigma_zs, strain #Gives the normal stress function for a given span zs, and x- and y- coordinate

# print('sigma_zs', Normal_stress_due_to_bending(0.15, Wing.airfoilordinate(Wing.c)))



def Shear_wb(zs, dL, dD, dM):
    #section 01
    qtorque = dM/(2*Wing.Area_cell())
    n = 100
    ds = Wing.HSpar1/n
    qs1 = np.array([])
    s1 = np.array([])
    qs1 = np.append(qs1, 0)
    s1 = np.append(s1, 0)
    s = Q_(ureg("0 m"))
    for i in range(n):
        s = s + ds
        s1 = np.append(s1, s)
        qs = s**2*Wing.ThSpar1*(-dL)/Inertia.Ixx_wb
        qs += s*Wing.ChSpar1*Wing.length_chord(zs)*Wing.ThSpar1*(-dD)/Inertia.Iyy_wb
        qs += qtorque
        qs1 = np.append(qs1, qs)
    section01at1 = qs1[-1]
    section01at1 *= Q_("N/m")

    #section12
    n = 100 #number of sections
    ds = (Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)/n)
    qs2 = np.array([])
    s2 = np.array([])
    qs2 = np.append(qs2, section01at1)
    s2 = np.append(s2, 0)
    s = 0
    for i in range(n):
        s = s + ds
        s2 = np.append(s2, s)
        qs = s * Wing.get_xy_from_perim(s/Wing.length_chord(zs))[0]*Wing.length_chord(zs)*Wing.ThSkin*(-dL)/Inertia.Ixx_wb
        qs += s * Wing.get_xy_from_perim(s/Wing.length_chord(zs))[1]*Wing.length_chord(zs)*Wing.ThSkin*(-dD)/Inertia.Iyy_wb
        qs +=  section01at1
        qs2 = np.append(qs2, qs)
    section12at2 =  qs2[-1]
    section12at2 *= Q_("N/m")
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
        qs = -s**2*Wing.ThSpar2*(-dL)/Inertia.Ixx_wb
        qs += -s*Wing.ChSpar2*Wing.length_chord(zs)*Wing.ThSpar1*(-dD)/Inertia.Iyy_wb
        qs += section12at2
        qs3 = np.append(qs3, qs)
    qs1 *= ureg("N/m")
    qs2 *= ureg("N/m")
    qs3 *= ureg("N/m")
    s1 *= ureg("m")
    s2 *= ureg("m")
    s3 *= ureg("m")
    #print("s1=,", s1)
    #print("q1=", qs1)
    #print("s2=,", s2)
    #print("q2=", qs2)
    #print("s3=,", s3)
    #print("q3=", qs3)
    return s1, s2, s3, qs1, qs2, qs3




def calc_qs0(Shear_wb, dL, zs, dD, dM):
    n = 101
    ds = Wing.HSpar1/n
    s1 = Shear_wb(zs, dL, dD, dM)[0]
    s2 = Shear_wb(zs, dL, dD, dM)[1]
    s3 = Shear_wb(zs, dL, dD, dM)[2]
    qs1 = Shear_wb(zs, dL, dD, dM)[3]
    qs2 = Shear_wb(zs, dL, dD, dM)[4]
    qs3 = Shear_wb(zs, dL, dD, dM)[5]
    qs0nom = 0
    for i in range(n):
        qs0nom += qs1[i]* ds
    ds = (Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)/n)
    for i in range(n):
        qs0nom += qs2[i]* ds
    ds = Wing.HSpar2/n
    for i in range(n):
        qs0nom += qs3[i]* ds
    qs0denom = Wing.HSpar1/Wing.ThSpar1 + Wing.HSpar2/Wing.ThSpar2 + 2*Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)/Wing.ThSkin
    qs0denom *= Q_(ureg("m"))
    qs0 =  - qs0nom/qs0denom
    for i in range(n):
        qs1[i] += qs0
        qs2[i] += qs0
        qs3[i] += qs0
    return s1, s2, s3, qs1, qs2, qs3

s1, s2, s3, qs1, qs2, qs3 = calc_qs0(Shear_wb, dL, zs, dD, dM)

#print("s1=,", s1)
print("q1=", qs1)
#print("s2=,", s2)
print("q2=", qs2)
#print("s3=,", s3)
print("q3=", qs3)

#Moment_from_shear = calc_moment_from_shear(Wing.z,s1, s2, s3, qs1, qs2, qs3)

def calc_moment_from_shear(zs, s_1, s_2, s_3, q_1, q_2, q_3):
    Moment = 0
    Momentz = np.array([])
    if (len(q_2) != len(s_2)):
        raise ValueError("ERROR, ARRAY LENGTH NOT EQUAL!")
    for i in range(0,len(s_2)-1,2):
        ds = s_2[i+1] - s_2[i]
        q_loc = (q_2[i] + q_2[i+1])*0.5
        x_loc_1, y_loc_1 = Wing.get_xy_from_perim(s_2[i]/Wing.length_chord(zs), Wing.ChSpar1)
        x_loc_1 += -1*Wing.ChSpar1
        x_loc_2, y_loc_2 = Wing.get_xy_from_perim(s_2[i+1]/Wing.length_chord(zs), Wing.ChSpar1)
        x_loc_2 += -1 * Wing.ChSpar1
        x_loc_1 *= Wing.length_chord(zs)
        y_loc_1 *= Wing.length_chord(zs)
        x_loc_2 *= Wing.length_chord(zs)
        y_loc_2 *= Wing.length_chord(zs)
        slope = (y_loc_2 - y_loc_1)/(x_loc_2 - x_loc_1)
        force_angle = np.arctan(slope)
        # print("q_loc=", q_loc)
        # print("ds=", ds)
        Fx = q_loc*ds*np.cos(force_angle)
        Fy = q_loc*ds*np.sin(force_angle)
        # print("Fx=", Fx)
        # print("Fy=", Fy)
        Moment += -Fx*(y_loc_1+y_loc_2)*0.5
        Moment += Fy*(x_loc_1+x_loc_2)*0.5
        Momentz = np.append(Momentz, Moment)

    for i in range(0, len(s_3)-2, 2):
        ds = s_3[i + 1] - s_3[i]
        q_loc = (q_3[i] + q_3[i + 1]) * 0.5
        x_loc = Wing.ChSpar2*Wing.length_chord(zs)
        # print("q_loc=", q_loc)
        # print("ds=", ds)
        # print("Fy=", q_loc * ds)
        Moment += -x_loc * q_loc * ds
        Momentz = np.append(Momentz, Moment)

    # plt.plot(Momentz)
    # plt.show()
    return Moment
Moment = calc_moment_from_shear(Wing.z,s1, s2, s3, qs1, qs2, qs3)

print("Moment from shear", Moment)

def Shearcenter(Moment, dL):
    shearcentre_x = Moment/dL
    return shearcentre_x

shearcentre_x = Shearcenter(Moment, dL)
print("shearcenter", shearcentre_x)


def Torsion(dM, dL, shearcentre_x):
    A_cell = Wing.Area_cell()
    T = dM + dL*shearcentre_x
    const_tor = T/(4*A_cell**2*shear_modulus) #constant term in twist formula
    print("0", const_tor)
    line_int_tor  = 2*Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)*Wing.length_chord(zs)/Wing.ThSkin
    line_int_tor = Wing.HSpar1*Wing.length_chord(zs)/Wing.ThSpar1
    line_int_tor = Wing.HSpar2*Wing.length_chord(zs)/Wing.ThSpar2 #result from line integral from torsion formula
    twist_wb_tor_per_m  = const_tor*line_int_tor
    return twist_wb_tor_per_m

print("twist", Torsion(dM, dL, shearcentre_x))

#Final shear flow
    
def Shear_stress_wb(zs, dL, dD, dM):
    n = 100
    s1 = calc_qs0[0]
    s2 = calc_qs0[1]
    s3 = calc_qs0[2]
    qs1 = calc_qs0[3]
    qs2 = calc_qs0[4]
    qs3 = calc_qs0[5]
    tau1 = np.array([])
    tau2 = np.array([])
    tau3 = np.array([])
    qs1 = np.append(qs1, 0)
    #section 01
    for i in range(n):
        tauinrange1 = qs1[i]/Wing.ThSpar1
        tau1 = np.append(tau1, tauinrange1) #tau in y direction
    #section12
    for i in range(n):
        tauinrange2 = qs2[i]/Wing.ThSkin
        tau2 = np.append(tau2, tauinrange2)
        
    #section23
    for i in range(n):
        tauinrange3 = qs3[i]/Wing.ThSpar2
        tau3 = np.append(tau3, tauinrange3) #tau in y direction
    tau1 *= ureg("N/m**2")
    tau2 *= ureg("N/m**2")
    tau3 *= ureg("N/m**2")
    s1 *= ureg("m")
    s2 *= ureg("m")
    s3 *= ureg("m")

    return s1, s2, s3, tau1, tau2, tau3

# print("twist =", Torsion(Geometry.Wing.b/2,calc_qbase(Moment_from_shear, L, Wing.z), M).to("rad/m"))

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
def Tsia_Wu(sigma_zs):
    F11=1/(yield_strength*compr_strength)
    F22 = F11
    F12 = -1/2*np.sqrt(F11*F22)
    F1 = 1/(yield_strength)-1/(compr_strength)
    F2 = 1/(yield_strength)-1/(compr_strength)
    F44 = 1/tau_max**2
    F66 = 1/tau_max**2
    sigma1 = sigma_zs
    sigma2 = 0
    sigma3 = 0
    tau12 = 1 #DUMMY VALUE
    tau23 = 0 #DUMMY VALUE
    tau13 = 1 #DUMMY VALUE 
    F = F11 *sigma1**2+F22*(sigma2**2+sigma3**2)+sigma2*sigma3*(2*F22-F44)
    F += 2*F12*sigma1*(sigma3+sigma2)+F1*(sigma1+sigma2) + F2*sigma3
    F += F44*tau23**2 + F66*(tau13**2+tau12**2)
    if F < 1:
        print("No failure occurs")
    else:
        print("Failure occurs")



#plt.plot(s3, qs3)
#plt.show()

# with is like your try .. finally block in this case
with open('StrucVal.py', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

data[18] = 'youngs_modulus = Q_(\"' + str(youngs_modulus) + '\")\n'

# now change the 2nd line, note that you have to add a newline

# and write everything back
with open('StrucVal.py', 'w') as file:
    file.writelines(data)

