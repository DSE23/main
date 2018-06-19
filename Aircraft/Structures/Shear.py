"""
Name: Shear
Department: Structures
Last updated: 19/06/2018 12:30 by Midas
"""

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_ # Imports the unit registry from the Misc folder
from Structures import WingStress
from Structures import Wing
from Structures import Inertia
import numpy as np

# Calculate Boom Area's [Midas]
# Units checked and correct
def Get_boom_area(A_spar_caps):
    """
    This function calculates the boom area's needed for the structural idealization
    It assumes the wingbox is loaded in bending.
    :param A_spar_caps: Spar cap area's [Expected With Unit]
    :return: Boom area array with dimension.
    """
    x_coords = Wing.x_y_angle_coords[0]
    y_coords = Wing.x_y_angle_coords[1]
    n_str = len(x_coords)
    t_sk = Wing.ThSkin
    h_sp_1 = Wing.HSpar1
    h_sp_2 = Wing.HSpar2
    B_area = np.array([])

    B_1 = A_spar_caps #+ (Wing.ThSpar1*h_sp_1/6)*(2 + (-1)) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[1]/y_coords[0])
    B_area = np.append(B_area, B_1.to(ureg("m**2")))

    B_2 = Wing.A_stringer #+ (t_sk*Wing.perim_spacing/6)*(2 + h_sp_1/2/ y_coords[0]) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[1]/ y_coords[0])
    B_area = np.append(B_area, B_2.to(ureg("m**2")))
    for i in range(1, n_str-1):
        B_i = Wing.A_stringer #+ (t_sk*Wing.perim_spacing/6)*(2 + y_coords[i-1]/ y_coords[i]) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[i+1]/ y_coords[i])
        B_area = np.append(B_area, B_i.to(ureg("m**2")))

    B_n_min_1 = Wing.A_stringer #+ (t_sk*Wing.perim_spacing/6)*(2 + h_sp_2/2 / y_coords[-1]) + (t_sk*Wing.perim_spacing/6)*(2 + y_coords[-2]/ y_coords[-1])
    B_area = np.append(B_area, B_n_min_1.to(ureg("m**2")))

    #B_n = (t_sk*Wing.perim_spacing/6)*(2 + y_coords[-2]/ y_coords[-1]) + (Wing.ThSpar2*h_sp_2/6)*(2 + (-1))
    #B_area = np.append(B_area, B_n.to(ureg("m**2")))

    B_area = np.append(B_area, np.flip(B_area, 0))

    return B_area*ureg("m**2")
# UNCOMMENT TO TEST GET BOOM AREA FUNCTION:
print(Get_boom_area(Q_("1000 mm**2")))

# Calculate base shear flow for every section of the wing box # Midas & Tobias
def Calc_base_shear_flow(boom_areas):
    strs_x_coords, strs_y_coords, _ = Wing.x_y_angle_coords
    strs_x_coords.ito(ureg("m"))
    strs_y_coords.ito(ureg("m"))
    strs_x_coords = np.append(strs_x_coords, np.flip(strs_x_coords, 0))
    strs_x_coords *= ureg('m')
    strs_y_coords = np.append(strs_y_coords, np.flip(strs_y_coords, 0))
    strs_y_coords *= ureg('m')
    n = 100 # Number of sections
    S_x = WingStress.D
    S_y = WingStress.L
    Ixx = Inertia.Ixx_wb
    Iyy = Inertia.Iyy_wb
    ## section 1 2
    ds = Wing.HSpar1/(2*n)
    qs12L = np.array([])
    qs12D = np.array([])
    qs12L = np.append(qs12L, 0)
    qs12D = np.append(qs12D, 0)
    s1 = np.array([])
    s1 = np.append(s1, 0)
    s = Q_("0 m")
    q_loc_L = 0
    q_loc_D = 0
    x = (Wing.ChSpar1 - Wing.centroid) * Wing.Chordlength  # x_coordinate of Spar 1 w.r.t. the centroid
    for _ in range(n-1):
        s += ds
        y = s
        s1 = np.append(s1, s)
        # Lift
        q_loc_L += (S_y/Ixx)*Wing.ThSpar1*y*ds
        qs12L = np.append(qs12L, q_loc_L)
        # Drag
        q_loc_D += (S_x/Iyy)*Wing.ThSpar1*x*ds
        qs12D = np.append(qs12D, q_loc_D)
    # Lift
    q_loc_L += (S_y / Ixx) *( Wing.ThSpar1 * y * ds + boom_areas[0] * Wing.HSpar1/2 )
    qs12L = np.append(qs12L, q_loc_L)
    # Drag
    q_loc_D += (S_x / Iyy) *( Wing.ThSpar1 * x * ds + boom_areas[0] * x)
    qs12D = np.append(qs12D, q_loc_D)

    qs12L *= ureg("N/m")
    qs12D *= ureg("N/m")
    ## section 2 3
    ds = (Wing.skin_length)/n
    s = Q_("0 m")
    qs23L = np.array([])
    qs23L = np.append(qs23L, qs12L[-1])
    qs23D = np.array([])
    qs23D = np.append(qs23D, qs12D[-1])
    s2 = np.array([])
    s2 = np.append(s2, s)
    str_counter = 0
    for _ in range(n):
        s = np.add(ds, s)
        s2 = np.append(s2, s)
        x_coor, y_coor = Wing.get_xy_from_perim(s/Wing.Chordlength, Wing.ChSpar1) # RELATIVE!!
        # Lift
        q_loc_L += -(S_y/Ixx)*(Wing.ThSkin*y_coor*Wing.Chordlength*ds)

        # Drag
        q_loc_D += -(S_x / Iyy) * (Wing.ThSkin * (x_coor - Wing.centroid) * Wing.Chordlength * ds)
        if(abs(x_coor*Wing.Chordlength - strs_x_coords[str_counter]) < Q_("10e-3 mm")):
            q_loc_L +=  Wing.A_stringer*strs_y_coords[str_counter]
            q_loc_D += Wing.A_stringer *(strs_x_coords[str_counter]/Wing.Chordlength - Wing.centroid)*Wing.Chordlength
            str_counter += 1
        qs23L = np.append(qs23L, q_loc_L)
        qs23D = np.append(qs23D, q_loc_D)

    ## Section 3 5
    ds = Wing.HSpar2/n
    qs35L = np.array([])
    qs35D = np.array([])
    qs35L = np.append(qs35L, qs23L[-1])
    qs35D = np.append(qs35D, qs23D[-1])
    s3 = 0
    s3 = np.append(s3, 0)
    s = Q_("0 m")
    x = (Wing.ChSpar2 - Wing.centroid) * Wing.Chordlength  # x_coordinate of Spar 1 w.r.t. the centroid

    for _ in range(n):
        s += ds
        y = Wing.HSpar2/2 - s
        s3 = np.append(s3, s)
        # Lift
        q_loc_L += (S_y / Ixx) * Wing.ThSpar2 * y * ds
        qs35L = np.append(qs35L, q_loc_L)
        # Drag
        q_loc_D += (S_x / Iyy) * Wing.ThSpar2 * x * ds
        qs35D = np.append(qs35D, q_loc_D)

    qs35L *= ureg("N/m")
    qs35D *= ureg("N/m")

    ## section 5 6
    ds = (Wing.skin_length) / n
    s = Q_("0 m")
    qs56L = np.array([])
    qs56L = np.append(qs56L, qs35L[-1])
    qs56D = np.array([])
    qs56D = np.append(qs56D, qs35D[-1])
    s4 = np.array([])
    s4 = np.append(s4, s)
    for _ in range(n):
        s = np.add(ds, s)
        s4 = np.append(s4, s)
        x_coor, y_coor = Wing.get_xy_from_perim(s / Wing.Chordlength, Wing.ChSpar2, reverse=True)  # RELATIVE!!
        # Lift
        q_loc_L += -(S_y / Ixx) * (Wing.ThSkin * y_coor * Wing.Chordlength * ds)

        # Drag
        q_loc_D += -(S_x / Iyy) * (Wing.ThSkin * (x_coor - Wing.centroid) * Wing.Chordlength * ds)
        if (abs(x_coor * Wing.Chordlength - strs_x_coords[str_counter]) < Q_("10e-3 mm")):
            q_loc_L += Wing.A_stringer * strs_y_coords[str_counter]
            q_loc_D += Wing.A_stringer * (
                        strs_x_coords[str_counter] / Wing.Chordlength - Wing.centroid) * Wing.Chordlength
            str_counter += 1
        qs56L = np.append(qs56L, q_loc_L)
        qs56D = np.append(qs56D, q_loc_D)

    ## section 6 1
    ds = Wing.HSpar1/(2*n)
    qs61L = np.array([])
    qs61D = np.array([])
    qs61L = np.append(qs61L, q_loc_L)
    qs61D = np.append(qs61D, q_loc_D)
    s5 = np.array([])
    s5 = np.append(s5, 0)
    s = Q_("0 m")
    x = (Wing.ChSpar1 - Wing.centroid) * Wing.Chordlength  # x_coordinate of Spar 1 w.r.t. the centroid
    # Lift
    q_loc_L += (S_y / Ixx) *( Wing.ThSpar1 * y * ds + boom_areas[-1] * Wing.HSpar1/2 )
    qs61L = np.append(qs12L, q_loc_L)
    # Drag
    q_loc_D += (S_x / Iyy) *( Wing.ThSpar1 * x * ds + boom_areas[-1] * x)
    qs61D = np.append(qs12D, q_loc_D)
    for _ in range(n-1):
        s += ds
        y = -Wing.HSpar1/2 + s  # x_coordinate of Spar 1 w.r.t. the centroid
        s5 = np.append(s5, s)
        # Lift
        q_loc_L += (S_y/Ixx)*Wing.ThSpar1*y*ds
        qs61L = np.append(qs61L, q_loc_L)
        # Drag
        q_loc_D += (S_x/Iyy)*Wing.ThSpar1*x*ds
        qs61D = np.append(qs61D, q_loc_D)


    qs12L *= ureg("N/m")
    qs12D *= ureg("N/m")

    return qs12L, qs23L, qs35L

print(Calc_base_shear_flow(Get_boom_area(Q_("1000 mm**2"))))



# Calculate correcting shear flow (qs0)      #Tobias #WORK IN PROGRESS
def calc_correcting_shear_flow():
#s1, s2, s3, s4, s5, qs1L, qs2L, qs3L, qs4L, qs5L, qs1D, qs2D, qs3D, qs4D, qs5D  = b
#    for i in range(5):

    #qs0nom =
    qs0denom = Wing.HSpar1/Wing.ThSpar1
    qs0denom += 2*Wing.length_Skin_x_c/Wing.ThSkin
    qs0denom += Wing.HSpar2/Wing.ThSpar2
    qs0 = -qs0nom/qs0denom
    return qs0

# Add correcting shear flow to base shear flows
    
# Add Moment shear flow to base shear flows



# Compute moments around a.c. caused by shear forces due to shear flows

# Calculate shear center location   #Tobias
#units checked and correct
def Shear_center(moment_shear):
    shear_center = moment_shear/WingStress.L
    return shear_center


# Calculate Torque                 #Tobias
#units checked and correct
def Torque_for_twist(shear_center):
    T = WingStress.M + WingStress.L * shear_center
    return T


# Calculate Rate of Twist          #Tobias
#units checked and correct
def Rate_of_twist(T):
    constant = T/(4*Wing.Area_cell()**2*WingStress.shear_modulus)
    integral = Wing.HSpar1/Wing.ThSpar1
    integral += 2*Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)/Wing.ThSkin
    integral += Wing.HSpar2/Wing.ThSpar2
    dthetadz = constant/integral
    return dthetadz

# Calculate shear stress