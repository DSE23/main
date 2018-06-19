"""
Name: Shear
Department: Structures
Last updated: 19/06/2018 10:04 by Midas
"""

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_ # Imports the unit registry from the Misc folder
from Structures import WingStress
from Structures import Wing
import numpy as np
# Calculate Boom Area's
def get_boom_area(A_spar_caps):
    x_coords = Wing.x_y_angle_coords[0]
    y_coords = Wing.x_y_angle_coords[1]
    n_str = len(x_coords)
    t_sk = Wing.ThSkin
    h_sp_1 = Wing.HSpar1
    h_sp_2 = Wing.HSpar2
    B_area = np.array([])

    B_1 = A_spar_caps + (Wing.ThSpar1*h_sp_1/6)*(2 + (-1)) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[1]/y_coords[0])
    B_area = np.append(B_area, B_1.to(ureg("m**2")))

    B_2 = Wing.A_stringer + (t_sk*Wing.perim_spacing/6)*(2 + h_sp_1/2/ y_coords[0]) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[1]/ y_coords[0])
    B_area = np.append(B_area, B_2.to(ureg("m**2")))
    for i in range(1, n_str-1):
        B_i = Wing.A_stringer + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[i-1]/ y_coords[i]) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[i+1]/ y_coords[i])
        B_area = np.append(B_area, B_i.to(ureg("m**2")))

    B_n_min_1 = Wing.A_stringer + (t_sk*Wing.perim_spacing/6)*(2 + h_sp_2/2 / y_coords[-1]) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[-2]/ y_coords[-1])
    B_area = np.append(B_area, B_n_min_1.to(ureg("m**2")))

    B_n = (t_sk*Wing.perim_spacing/6)*(2 + y_coords[-2]/ y_coords[-1]) + (Wing.ThSpar2*h_sp_2/6)*(2 + (-1))
    B_area = np.append(B_area, B_n.to(ureg("m**2")))

    B_area = np.append(B_area, np.flip(B_area, 0))

    return B_area*ureg("m**2")
# UNCOMMENT TO TEST GET BOOM AREA FUNCTION:
#print(get_boom_area(Q_("1000 mm**2")))

# Calculate base shear flow for every section of the wing box

# Calculate correcting shear flow (qs0)      #Tobias
def calc_corr_shearflow():
#s1, s2, s3, s4, s5, qs1L, qs2L, qs3L, qs4L, qs5L, qs1D, qs2D, qs3D, qs4D, qs5D  = b
#    for i in range(5):

    #qs0nom =
    qs0denom = Wing.HSpar1/Wing.ThSpar1
    qs0denom += 2*Wing.length_Skin_x_c/Wing.ThSkin
    qs0denom += Wing.HSpar2/Wing.ThSpar2
    qs0 = -qs0nom/qs0denom
    return qs0

# Add correcting shear flow to base shear flows

# Compute moments around a.c. caused by shear forces due to shear flows

# Calculate shear center location   #Tobias
def Shear_center(moment_shear):
    shear_center = moment_shear/WingStress.L
    return shear_center

print(Shear_center(Q_("200000 N*m")))

# Calculate Torque                 #Tobias
def Torque_for_twist(shear_center):
    T = WingStress.M + WingStress.L * shear_center
    return T

# Calculate Rate of Twist          #Tobias
def Rate_of_twist(T):
    constant = T/(4*Wing.Area_cell()*WingStress.shear_modulus)
    integral = Wing.HSpar1/Wing.ThSpar1
    integral += 2*Wing.length_Skin_x_c/Wing.ThSkin
    integral += Wing.HSpar2/Wing.ThSpar2
    dthetadz = constant/integral
    return dthetadz

# Calculate shear stress