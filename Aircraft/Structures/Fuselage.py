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
rho = Performance.rho_c.magnitude * ureg("kg/(m**3)") #rho at cruise altitude
V = Performance.V_cruise.magnitude * ureg("m/s")        #cruise speed
S_VT = Geometry.V_Tail.S                    #area of vertical tail
S_HT = Geometry.H_Tail.S                    #are of horizontal tail
b_VT = Geometry.V_Tail.b                    #span vertical tail
MAC_VT = Geometry.V_Tail.MAC
MAC_HT = Geometry.H_Tail.MAC

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

'''------------------Inertia calculation-----------------------'''

#Inertia section 1 (normal)
Izz_sec1, Iyy_sec1 = (b_f80/2)**2 * B_sec1 * 4

#Inertia section 2
Izz_sec2, Iyy_sec2 = (b_f80/2)**2 * B_sec2 * 4

#Inertia section 3
Izz_sec3, Iyy_sec3 = (b_f_taper/2)**2 * B_sec3 * 4

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
Mx_sec2 = (b_VT * 0.33) * L_VT
Mz_sec2 = (l_fus - x) * L_HT - (b_VT * 0.33) * D_VT

#Bending Loading section 3
My_sec3 = (l_fus - x) * L_VT
Mx_sec3 = (b_VT * 0.33) * L_VT
Mz_sec3 = (l_fus - x) * L_HT - (b_VT * 0.33) * D_VT


#Axial Loading of section 1, 2 & 3
ax_sec1 = Engine_mounts.f_x
ax_sec2 = D_VT + D_HT
ax_sec3 = D_VT + D_HT


'''----------Bending stress calculations----------------'''

'''section 1'''
#Bending section 1
sigma_x_sec1 = My_sec1 / Iyy_sec1 * z + Mz_sec1 / Izz_sec1 * y + (ax_sec1 / (B_sec1 * 4))        #bending stress

#Torsion section 1

q_x_sec1 = Mx_sec1 / (2 * b_f80 * 2)                                  #shear flow
shear_x_sec1 = q_x_sec1 / t

'''section 2'''
#Bending section 2
sigma_x_sec2 = My_sec2 / Iyy_sec2 * z + Mz_sec2 / Izz_sec2 * y + (ax_sec2 / (B_sec2 * 4))

#Torsion section 2

q_x_sec2 = Mx_sec2 / (2 * b_f80 * 2)
shear_x_sec2 = q_x_sec2 / t

'''section 3'''
#Bending section 3
sigma_x_sec3 = My_sec3 / Iyy_sec3 * z + Mz_sec3 / Izz_sec3 * y + (ax_sec3 / (B_sec3 * 4))

#Torsion section 3
q_x_sec3 = Mx_sec3 / (2 * b_f_taper * 2)
shear_x_sec3 = q_x_sec3 / t



print(sigma_x_sec1)










