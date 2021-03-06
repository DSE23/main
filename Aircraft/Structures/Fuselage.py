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
from Propulsion_and_systems import Engine_mounts
from Aerodynamics import Wing as AWing
from Performance import Performance
import matplotlib
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D

from Misc import ureg, Q_ # Imports the unit registry from the Misc folder

b_f = Geometry.Fuselage.b_f                   #diameter of fuselage at the fire wall
b_f80 = b_f*0.8                               #design radius
print('b_f80', b_f80)
t = Q_('0.0025 m')                             #thickness of fuselage skin
t_ribs = Q_('0.002 m')                          #thickness of the ribs
h_ribs = Q_('0.05 m')                           #the height of the ribs
s_ribs = Q_('0.5 m')                            #rib spacing
moter_w = Q_('180 kg')                          #Resultant weight of the motor
g = Q_('9.81 m / s**2')                         #The gravity acceleration
l_fus = Geometry.Fuselage.l_f                   #length of the fuselage in m
print('l_fus', l_fus)
l_sec1 = Q_('1.5 m')                          #length of section 1 (normal)
l_sec2 = Q_('1.0 m')                          #length of section 2 (cut out)
l_sec3 = l_fus - l_sec1 - l_sec2               #length of section 3 (taper)
print('l_sec3', l_sec3)
rho = Performance.rho_c.magnitude * ureg("kg/(m**3)") #rho at cruise altitude
V = Performance.V_cruise.magnitude * ureg("m/s")        #cruise speed
S_VT = Geometry.V_tail.S                    #area of vertical tail
S_HT = Geometry.H_tail.S                    #are of horizontal tail
b_VT = Geometry.V_tail.b                    #span vertical tail
MAC_VT = Geometry.V_tail.MAC
MAC_HT = Geometry.H_tail.MAC

#Material properties of the chosen material.
#Current chosen material:
#Epoxy/Carbon fiber, UD prepreg, QI lay-up
youngs_modulus = Q_("60.1 GPa")  #E
yield_strength = Q_("738 MPa")  #tensile
compr_strength = Q_("657 MPa") #compression
shear_modulus = Q_("23 GPa")   #G
poisson = 0.31                 # maximum 0.33
tau_max = Q_("40 MPa")
density = Q_("1560 kg/m**3")

'''Not very elegant solution for problem'''
# Boom areas section 1 (normal)
B_sec1 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2  # area for all booms

# Boom areas section 2 (cut out)
B_sec2 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2


def B_calc(x):
    if x <= l_sec1 + l_sec2:
        b_f_taper = b_f80
        B_sec3 = B_sec1

    if x > l_sec1 + l_sec2:
        b_f_taper = b_f80 - ((b_f80 / (l_sec3 + l_sec2)) * (x - l_sec2 - l_sec1))
        B_sec3 = t * (b_f_taper / 2) / 2 + t * (b_f_taper / 2) / 2

    return B_sec3, b_f_taper



'''The definition that caclulates certain normal and shear stresses at specific points in the fuselage'''

def fuselage_calc(x, y, z):

    #Boom areas section 1 (normal)
    B_sec1 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2          #area for all booms


    #Boom areas section 2 (cut out)
    B_sec2 = t * (b_f80 / 2) / 2 + t * (b_f80 / 2) / 2

    #Boom areas section 3 (cone/taper)
    def B_calc(x):
        if Q_('0 m') <= x <= l_sec1 + l_sec2:
            b_f_taper = b_f80
            B_sec3 = B_sec1

        if l_sec1 + l_sec2 < x <= l_fus:
            b_f_taper = b_f80 - ((b_f80 / (l_sec3 + l_sec2)) * (x - l_sec2 - l_sec1))
            B_sec3 = t * (b_f_taper / 2) / 2 + t * (b_f_taper / 2) / 2


        return B_sec3, b_f_taper

    # y = B_calc(x)[1]
    # z = B_calc(x)[1]                                         #the dreadful z


    B_sec3, b_f_taper = B_calc(x)

    '''------------------Inertia calculation-----------------------'''

    #Inertia section 1 (normal)
    Iyy_sec1 = (b_f80/2)**2 * B_sec1 * 4
    Izz_sec1 = Iyy_sec1

    #Inertia section 2
    Iyy_sec2 = (b_f80/2)**2 * B_sec2 * 4
    Izz_sec2 = Iyy_sec2

    #Inertia section 3
    Iyy_sec3 = (b_f_taper/2)**2 * B_sec3 * 4
    Izz_sec3 = Iyy_sec3

    print(x, b_f80, b_f_taper)
    '''------------Loading force/moment calculation------------'''

    #Tail force calculation
    cl_vt, cd_vt, cm_vt = AWing.computeloadsvt()
    cl_ht, cd_ht, cm_ht = AWing.computeloadsht()
    L_VT = 0.5 * cl_vt * rho * (V ** 2) * S_VT                          #Tangent forces of the control surfaces at full deflection
    L_HT = 0.5 * cl_ht * rho * (V ** 2) * S_HT

    D_VT = 0.5 * cd_vt * rho * (V ** 2) * S_VT                          #Normal forces of the control surfaces at full deflection
    D_HT = 0.5 * cd_ht * rho * (V ** 2) * S_HT

    M_VT = 0.5 * cm_vt * rho * (V ** 2) * S_VT * MAC_VT
    M_HT = 0.5 * cm_vt * rho * (V ** 2) * S_HT * MAC_HT


    #Bending Loading section 1
    My_sec1 = Engine_mounts.m_y         #calculation of the resultant moment in y as variable of z
    Mx_sec1 = Engine_mounts.m_x         #calculation of the resultant moment in y as variable of z
    Mz_sec1 = x * Engine_mounts.f_z + Engine_mounts.m_z          #calculation of the resultant moment in y as variable of z

    #Bending Loading section 2
    My_sec2 = (l_fus - x) * L_VT
    My_sec2.ito(ureg('N * m'))
    Mx_sec2 = (b_VT * 0.33) * L_VT
    Mx_sec2.ito(ureg('N * m'))
    Mz_sec2 = -(l_fus - x) * L_HT - (b_VT * 0.33) * D_VT
    Mz_sec2.ito(ureg('N * m'))

    #Bending Loading section 3
    My_sec3 = My_sec2
    Mx_sec3 = Mx_sec2
    Mz_sec3 = Mz_sec2


    #Axial Loading of section 1, 2 & 3
    ax_sec1 = Engine_mounts.f_x
    ax_sec2 = D_VT + D_HT
    ax_sec3 = ax_sec2


    '''----------Bending stress calculations----------------'''



    def normal_shear_stress(x):

        if Q_('0 m') <= x <= l_sec1:
            '''section 1'''

            #Bending section 1
            sigma_x = My_sec1 / Iyy_sec1 * z + Mz_sec1 / Izz_sec1 * y + (ax_sec1 / (B_sec1 * 4))        #bending stress
            #Torsion section 1

            q_x_sec1 = Mx_sec1 / (2 * b_f80 * b_f80)                                  #shear flow
            shear_x = q_x_sec1

        if l_sec1 < x < l_sec1 + l_sec2:

            '''section 2'''
            #Bending section 2
            sigma_x = My_sec2 / Iyy_sec2 * z + Mz_sec2 / Izz_sec2 * y + (ax_sec2 / (B_sec2 * 4))

            #Torsion section 2

            q_x_sec2 = Mx_sec2 / (2 * b_f80 * b_f80)
            shear_x = q_x_sec2

        if l_sec1 + l_sec2 <= x <= l_fus:
            '''section 3'''
            #Bending section 3
            sigma_x = My_sec3 / Iyy_sec3 * z + Mz_sec3 / Izz_sec3 * y + (ax_sec3 / (B_sec3 * 4))

            #Torsion section 3
            q_x_sec3 = Mx_sec3 / (2 * b_f_taper * b_f_taper)
            shear_x = q_x_sec3

        return sigma_x, shear_x

    sigma_x, shear_x = normal_shear_stress(x)

    '''------------Cut out correction calculation---------------'''
    q_12 = normal_shear_stress(x)[1]

    q_34_cor = (b_f80 * q_12) / b_f80
    q_23_cor = (b_f80 * q_34_cor) / b_f80
    q_14_cor = (b_f80 * q_23_cor) / b_f80
    print(q_34_cor, q_23_cor, q_14_cor)
    q_34 = q_12 + -q_34_cor
    q_23 = q_12 + q_23_cor
    q_14 = q_12 + q_14_cor
    q_12 = q_12 + -q_12


    P = ((q_14_cor + q_23_cor) * l_sec2) / 2


    P_stress = P / B_sec3

    if y > Q_('0 m') and z > Q_('0 m'):
        P_stress = -P_stress
    if y < Q_('0 m') and z < Q_('0 m'):
        P_stress = -P_stress

    if l_sec1 < x < l_sec1 + l_sec2:
        sigma_x = sigma_x + P_stress
        print('P_stress', P_stress)

        shear_x = q_34

        if z >= b_f_taper - Q_('0.01 m') or z <= -b_f_taper + Q_('0.01 m'):
            shear_x = q_23



    '''-------------Cut Out correction for non cut out parts (sec 1 and sec 3---------------'''

    if x <= l_sec1:
        a = np.array([[b_f_taper.magnitude, 0, -b_f_taper.magnitude, 0],
                     [0, b_f_taper.magnitude, 0, -b_f_taper.magnitude],
                     [0, 0, l_sec1.magnitude, l_sec1.magnitude],
                     [b_f_taper.magnitude, -b_f_taper.magnitude, 0, 0]])

        b = np.array([0, 0, P.magnitude, 0])

        M = np.linalg.solve(a, b)

        q_34 = shear_x + -M[2] * Q_('N/m')
        q_23 = shear_x + M[1] * Q_('N/m')
        q_14 = shear_x + M[3] * Q_('N/m')
        q_12 = shear_x + -M[0] * Q_('N/m')

        shear_x = q_12

        if z >= b_f_taper - Q_('0.01 m') or z <= -b_f_taper + Q_('0.01 m'):
            shear_x = q_23

    if x >= l_sec1 + l_sec2:
        a = np.array([[b_f_taper.magnitude, 0, -b_f_taper.magnitude, 0],
                      [0, b_f_taper.magnitude, 0, -b_f_taper.magnitude],
                      [0, 0, l_sec3.magnitude, l_sec3.magnitude],
                      [b_f_taper.magnitude, -b_f_taper.magnitude, 0, 0]])

        b = np.array([0, 0, P.magnitude, 0])

        M = np.linalg.solve(a, b)

        q_34 = shear_x + -M[2] * Q_('N/m')
        q_23 = shear_x + M[1] * Q_('N/m')
        q_14 = shear_x + M[3] * Q_('N/m')
        q_12 = shear_x + -M[0] * Q_('N/m')

        shear_x = q_12

        if z >= b_f_taper - Q_('0.01 m') or z <= -b_f_taper + Q_('0.01 m'):
            shear_x = q_23

    '''Area calculation'''
    Area = B_calc(x)[0]*4
    shear_x = shear_x / t

    '''Fuselage Ribs'''

    Area_ribs = b_f_taper * b_f_taper * h_ribs


    return sigma_x, shear_x, Area, Area_ribs


'''-----------------Tsia-Wu Failure criterion--------------------------'''
## For section 1
def Tsai_Wu(sigma_x, shear_x):
    F11=1/(yield_strength*compr_strength)
    F22 = F11
    F12 = -1/2*np.sqrt(F11*F22)
    F1 = 1/(yield_strength)-1/(compr_strength)
    F2 = 1/(yield_strength)-1/(compr_strength)
    F44 = 1/tau_max**2
    F66 = 1/tau_max**2
    sigma1 = sigma_x
    sigma2 = Q_("0 MPa")
    sigma3 = Q_("0 MPa")
    tau12 = shear_x
    tau23 = Q_("0 MPa")
    tau13 = Q_("0 MPa")
    F = F11 *sigma1**2+F22*(sigma2**2+sigma3**2)+sigma2*sigma3*(2*F22-F44)
    F += 2*F12*sigma1*(sigma3+sigma2)+F1*(sigma1+sigma2) + F2*sigma3
    F += F44*tau23**2 + F66*(tau13**2+tau12**2)
    F.ito(ureg('dimensionless'))
    if F.magnitude < 1:
        print("No failure occurs")
    if F.magnitude >= 1:
        print("Failure occurs")

    return F


'''-----------------------------Loop calcualtions----------------------------------'''

x = 0
n = 20                  #number of sections
ns = 4
x *= Q_('m')
sigmalist = np.array([])
shearlist = np.array([])
Flist = np.array([])
x_pos = np.array([])
y_pos = np.array([])
z_pos = np.array([])
Vol = Q_('0 m**3')
y = -B_calc(x)[1]
z = -B_calc(x)[1]
print('--------------', l_fus)

while x <= l_fus:
    print('next loop')
    print('x = ', x, 'l_fus =', l_fus)
    print('y = ', y)
    print('B_calc = ', B_calc(x)[1])
    y = - B_calc(x)[1]
    z = - B_calc(x)[1]
    while z <= B_calc(x)[1] - (B_calc(x)[1]/ns):
        print('1 x, y, z:', x, y, z)
        sigma_x, shear_x, Area, Area_ribs = fuselage_calc(x, y, z)
        F = Tsai_Wu(sigma_x, shear_x)
        sigmalist = np.append(sigmalist, sigma_x)
        shearlist = np.append(shearlist, shear_x)
        print('F', F)
        Flist = np.append(Flist, F.magnitude)
        x_pos = np.append(x_pos, x)
        y_pos = np.append(y_pos, y)
        z_pos = np.append(z_pos, z)
        # Vol = Vol + Area * (l_fus/n)
        print(sigma_x, shear_x)
        print('x, y, z', x, y, z)
        # Vol_rib = Area_ribs * t_ribs
        # N_ribs = float(int(l_fus / s_ribs))
        # Vol_ribs = Vol_rib * N_ribs
        z = z + (B_calc(x)[1] / ns)
    while y <= B_calc(x)[1] - (B_calc(x)[1]/ns):
        print('2 x, y, z:', x, y, z)
        sigma_x, shear_x, Area, Area_ribs = fuselage_calc(x, y, z)
        F = Tsai_Wu(sigma_x, shear_x)
        if l_sec1 < x < l_sec1 + l_sec2:
            if -B_calc(x)[1]+Q_('0.05 m') < y < B_calc(x)[1] - Q_('0.05 m'):
                sigma_x = Q_('0 Pa')
                F = Q_('0 dimensionless')
        sigmalist = np.append(sigmalist, sigma_x)
        shearlist = np.append(shearlist, shear_x)
        print('F', F)
        Flist = np.append(Flist, F.magnitude)
        x_pos = np.append(x_pos, x)
        y_pos = np.append(y_pos, y)
        z_pos = np.append(z_pos, z)
        # Vol = Vol + Area * (l_fus/n)
        print(sigma_x, shear_x)
        print('x, y, z', x, y, z)
        # Vol_rib = Area_ribs * t_ribs
        # N_ribs = float(int(l_fus / s_ribs))
        # Vol_ribs = Vol_rib * N_ribs
        y = y + (B_calc(x)[1] / ns)
    while z >= -B_calc(x)[1] + (B_calc(x)[1]/ns):
        print('3 x, y, z:', x, y, z)
        sigma_x, shear_x, Area, Area_ribs = fuselage_calc(x, y, z)
        F = Tsai_Wu(sigma_x, shear_x)
        sigmalist = np.append(sigmalist, sigma_x)
        shearlist = np.append(shearlist, shear_x)
        print('F', F)
        Flist = np.append(Flist, F.magnitude)
        x_pos = np.append(x_pos, x)
        y_pos = np.append(y_pos, y)
        z_pos = np.append(z_pos, z)
        # Vol = Vol + Area * (l_fus/n)
        print(sigma_x, shear_x)
        print('x, y, z', x, y, z)
        # Vol_rib = Area_ribs * t_ribs
        # N_ribs = float(int(l_fus / s_ribs))
        # Vol_ribs = Vol_rib * N_ribs
        z = z - (B_calc(x)[1] / ns)
    while y >= -B_calc(x)[1] + (B_calc(x)[1]/ns):
        print('4 x,y,z:', x, y, z)
        sigma_x, shear_x, Area, Area_ribs = fuselage_calc(x, y, z)
        F = Tsai_Wu(sigma_x, shear_x)
        sigmalist = np.append(sigmalist, sigma_x)
        shearlist = np.append(shearlist, shear_x)
        print('F', F)
        Flist = np.append(Flist, F.magnitude)
        x_pos = np.append(x_pos, x)
        y_pos = np.append(y_pos, y)
        z_pos = np.append(z_pos, z)
        # Vol = Vol + Area * (l_fus/n)
        print(sigma_x, shear_x)
        print('x, y, z', x, y, z)
        # Vol_rib = Area_ribs * t_ribs
        # N_ribs = float(int(l_fus / s_ribs))
        # Vol_ribs = Vol_rib * N_ribs
        y = y - (B_calc(x)[1] / ns)

    # while -B_calc(x - (l_fus/n))[1] <= y <= B_calc(x - +(l_fus/n))[1]:
    #     print('next loop')
    #     while -B_calc(x-(l_fus/n))[1] <= z <= B_calc(x-(l_fus/n))[1]:
    #         sigma_x, shear_x, Area, Area_ribs = fuselage_calc(x, y, z)
    #         F = Tsai_Wu(sigma_x, shear_x)
    #         sigmalist = np.append(sigmalist, sigma_x)
    #         shearlist = np.append(shearlist, shear_x)
    #         print('F', F)
    #         Flist = np.append(Flist, F.magnitude)
    #         x_pos = np.append(x_pos, x)
    #         y_pos = np.append(y_pos, y)
    #         z_pos = np.append(z_pos, z)
    #         # Vol = Vol + Area * (l_fus/n)
    #         print(sigma_x, shear_x)
    #         print('x, y, z', x, y, z)
    #         # Vol_rib = Area_ribs * t_ribs
    #         # N_ribs = float(int(l_fus / s_ribs))
    #         # Vol_ribs = Vol_rib * N_ribs
    #         z = z + (B_calc(x)[1] / 1)
    #     z = -B_calc(x)[1]
    #     y = y + (B_calc(x)[1] / 1)
    # y = -B_calc(x+(l_fus/n))[1]
    x = x + (l_fus / n)

Mass = Vol * density

print('Mass = ', Mass)

F = Flist
x = x_pos
y = y_pos
z = z_pos
s = sigmalist

cm = plt.get_cmap('rainbow')
if min(F) >= max(F):
    cNorm = matplotlib.colors.Normalize(vmin=min(F), vmax=max(F))
if min(F) < max(F):
    cNorm = matplotlib.colors.Normalize(vmin=min(F), vmax=max(F))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, c=scalarMap.to_rgba(F))
scalarMap.set_array(F)
fig.colorbar(scalarMap)
plt.show()

cm = plt.get_cmap('rainbow')
if min(s) >= max(s):
    cNorm = matplotlib.colors.Normalize(vmin=min(s), vmax=max(s))
if min(s) < max(s):
    cNorm = matplotlib.colors.Normalize(vmin=min(s), vmax=max(s))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, c=scalarMap.to_rgba(s))
scalarMap.set_array(s)
fig.colorbar(scalarMap)
plt.show()

plt.plot(x_pos, Flist)
plt.show()