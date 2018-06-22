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
from matplotlib import pyplot as plt



n = 400

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
#print(Get_boom_area(Q_("1000 mm**2")))


# Calculate base shear flow for every section of the wing box # Midas & Tobias
def Calc_base_shear_flow(boom_areas, n):
    """
    :param boom_areas: Return variable from  Boom_area function
    :param n: Number of sections that divides the perimiter
    :return: Arrays with s's and qs's, seperated for L and D
    """
    strs_x_coords, strs_y_coords, _ = Wing.x_y_angle_coords
    strs_x_coords.ito(ureg("m"))
    strs_y_coords.ito(ureg("m"))
    strs_x_coords = np.append(strs_x_coords, np.flip(strs_x_coords, 0))
    strs_x_coords *= ureg('m')
    strs_y_coords = np.append(strs_y_coords, np.flip(-1*strs_y_coords, 0))
    strs_y_coords *= ureg('m')
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
        q_loc_L += -(S_y/Ixx)*Wing.ThSpar1*y*ds
        qs12L = np.append(qs12L, q_loc_L.to(ureg("N/m")))
        # Drag
        q_loc_D += -(S_x/Iyy)*Wing.ThSpar1*x*ds
        qs12D = np.append(qs12D, q_loc_D.to(ureg("N/m")))
    s += ds
    y = s
    s1 = np.append(s1, s)
    # Lift
    q_loc_L += -(S_y / Ixx) *( Wing.ThSpar1 * y * ds + boom_areas[0] * Wing.HSpar1/2 )
    qs12L = np.append(qs12L, q_loc_L.to(ureg("N/m")))
    # Drag
    q_loc_D += -(S_x / Iyy) *( Wing.ThSpar1 * x * ds + boom_areas[0] * x)
    qs12D = np.append(qs12D, q_loc_D.to(ureg("N/m")))


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
        if(np.sqrt((x_coor*Wing.Chordlength - strs_x_coords[str_counter])**2 + (y_coor*Wing.Chordlength - strs_y_coords[str_counter])**2) < Q_("1 cm")):
            print("s2:", strs_x_coords[str_counter])
            q_loc_L += -(S_y/Ixx)*Wing.A_stringer*strs_y_coords[str_counter]
            q_loc_D += -(S_x/Iyy) *Wing.A_stringer *(strs_x_coords[str_counter]/Wing.Chordlength - Wing.centroid)*Wing.Chordlength
            str_counter += 1
        qs23L = np.append(qs23L, q_loc_L.to(ureg("N/m")))
        qs23D = np.append(qs23D, q_loc_D.to(ureg("N/m")))

    ## Section 3 5
    ds = Wing.HSpar2/n
    qs35L = np.array([])
    qs35D = np.array([])
    qs35L = np.append(qs35L, qs23L[-1])
    qs35D = np.append(qs35D, qs23D[-1])
    s3 = np.array([])
    s3 = np.append(s3, 0)
    s = Q_("0 m")
    x = (Wing.ChSpar2 - Wing.centroid) * Wing.Chordlength  # x_coordinate of Spar 1 w.r.t. the centroid

    for _ in range(n):
        s += ds
        y = Wing.HSpar2/2 - s
        s3 = np.append(s3, s)
        # Lift
        q_loc_L += -(S_y / Ixx) * Wing.ThSpar2 * y * ds
        qs35L = np.append(qs35L, q_loc_L.to(ureg("N/m")))
        # Drag
        q_loc_D += -(S_x / Iyy) * Wing.ThSpar2 * x * ds
        qs35D = np.append(qs35D, q_loc_D.to(ureg("N/m")))


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
        if(str_counter < 8):
            if (abs(x_coor * Wing.Chordlength - strs_x_coords[str_counter]) < Q_("1 cm")):
                q_loc_L += -(S_y / Ixx) * Wing.A_stringer * strs_y_coords[str_counter]
                print(strs_x_coords[str_counter])
                q_loc_D += -(S_x / Iyy) * Wing.A_stringer * (strs_x_coords[str_counter] / Wing.Chordlength - Wing.centroid) * Wing.Chordlength
                str_counter += 1
        qs56L = np.append(qs56L, q_loc_L.to(ureg("N/m")))
        qs56D = np.append(qs56D, q_loc_D.to(ureg("N/m")))

    ## section 6 1
    ds = Wing.HSpar1/(2*n)
    qs61L = np.array([])
    qs61D = np.array([])
    qs61L = np.append(qs61L, qs56L[-1])
    qs61D = np.append(qs61D, qs56D[-1])
    s5 = np.array([])
    s5 = np.append(s5, 0)
    s = Q_("0 m")
    x = (Wing.ChSpar1 - Wing.centroid) * Wing.Chordlength  # x_coordinate of Spar 1 w.r.t. the centroid
    y = -Wing.HSpar1 / 2
    # Lift
    s5 = np.append(s5, 0)
    q_loc_L += -(S_y / Ixx) *( Wing.ThSpar1 * y * ds + boom_areas[-1] * y )
    qs61L = np.append(qs61L, q_loc_L.to(ureg("N/m")))
    # Drag
    q_loc_D += -(S_x / Iyy) *( Wing.ThSpar1 * x * ds + boom_areas[-1] * x)
    qs61D = np.append(qs61D, q_loc_D.to(ureg("N/m")))
    for _ in range(n-1):
        s += ds
        y = -Wing.HSpar1/2 + s  # x_coordinate of Spar 1 w.r.t. the centroid
        s5 = np.append(s5, s)
        # Lift
        q_loc_L += -(S_y/Ixx)*Wing.ThSpar1*y*ds
        qs61L = np.append(qs61L, q_loc_L.to(ureg("N/m")))
        # Drag
        q_loc_D += -(S_x/Iyy)*Wing.ThSpar1*x*ds
        qs61D = np.append(qs61D, q_loc_D.to(ureg("N/m")))

    s1 *= ureg("m")
    s2 *= ureg("m")
    s3 *= ureg("m")
    s4 *= ureg("m")
    s5 *= ureg("m")
    qs12L *= ureg("N/m")
    qs12D *= ureg("N/m")
    qs23L *= ureg("N/m")
    qs23D *= ureg("N/m")
    qs35L *= ureg("N/m")
    qs35D *= ureg("N/m")
    qs56L *= ureg("N/m")
    qs56D *= ureg("N/m")
    qs61L *= ureg("N/m")
    qs61D *= ureg("N/m")

    return s1, s2, s3, s4, s5, qs12L, qs23L, qs35L, qs56L, qs61L, qs12D, qs23D, qs35D, qs56D, qs61D


s1, s2, s3, s4, s5, qs12L, qs23L, qs35L, qs56L, qs61L, qs12D, qs23D, qs35D, qs56D, qs61D  = Calc_base_shear_flow(Get_boom_area(Wing.AreaClamps/2), n)

#print("qs at 1st spar cap:", qs23L[0])
#print("qs at 2nd spar cap:", qs56L[-1])
#print("qs at beginning:", qs12L[0], "\tqs at end end:", qs61L[-1])

def Calculate_correcting_shear_flow(n):      #Tobias
    qs0denom = Wing.HSpar1/Wing.ThSpar1
    qs0denom += 2*Wing.skin_length/Wing.ThSkin
    qs0denom += Wing.HSpar2/Wing.ThSpar2
    qs0nomL = 0
    qs0nomD = 0
    for i in range(n+1):
        ds = s1[1]-s1[0]
        qs0nomL += qs12L[i]/Wing.ThSpar1*ds
        qs0nomD += qs12D[i]/Wing.ThSpar1*ds
        ds = s2[1]-s2[0]
        qs0nomL += qs23L[i]/Wing.ThSkin*ds
        qs0nomD += qs23D[i]/Wing.ThSkin*ds
        ds = s3[1]-s3[0]
        qs0nomL += qs35L[i]/Wing.ThSpar2*ds
        qs0nomD += qs35D[i]/Wing.ThSpar2*ds
        ds = s4[1]-s4[0]
        qs0nomL += qs56L[i]/Wing.ThSkin*ds
        qs0nomD += qs56D[i]/Wing.ThSkin*ds
        ds = s5[1]-s5[0]
        qs0nomL += qs61L[i]/Wing.ThSpar1*ds
        qs0nomD += qs61D[i]/Wing.ThSpar1*ds
    qs0L = -qs0nomL/qs0denom
    qs0D = -qs0nomD/qs0denom
    return qs0L, qs0D

qs0L, qs0D = Calculate_correcting_shear_flow(n)

# Add correcting shear flow to base shear flows
def Correcting_shearflow_array(n, qs0L, qs0D):
    qs0_L = np.array([])
    qs0_D = np.array([])
    for _ in range(n+1):
        qs0_L = np.append(qs0_L, qs0L)
        qs0_D = np.append(qs0_D, qs0D)
    return qs0_L*ureg("N/m"), qs0_D*ureg("N/m")
    
qs0_L, qs0_D = Correcting_shearflow_array(n, qs0L, qs0D)
# Add Moment shear flow to base shear flows

def Moment_shearflow(n):
    qmoment = WingStress.M/(2*Wing.Area_cell())
    q_moment = np.array([])
    for _ in range(n+1):
        q_moment = np.append(q_moment, qmoment)
    return q_moment*ureg("N/m")

q_moment = Moment_shearflow(n)
# Compute moments around a.c. caused by shear forces due to shear flows
def Calc_moment_due_to_shear(s1, s2, s3, s4, s5, qs12L, qs23L, qs35L, qs56L, qs61L, qs12D, qs23D, qs35D, qs56D, qs61D):

    t_ys12 = np.array([])

    x_coor_AC = 0.25*Wing.Chordlength
    Moment_L = 0
    # Moments from section 1 -> 2
    for i in range(0,len(qs12L)-1):
        q_loc = (qs12L[i] + qs12L[i+1])/2
        s_loc = (s1[i] + s1[i+1])/2
        ds = s1[i+1]-s1[i]
        x_loc = Wing.ChSpar1*Wing.Chordlength - x_coor_AC
        t_y = qs12L[i]/Wing.ThSpar1 + qs12D[i]/Wing.ThSpar1

        t_ys12 = np.append(t_ys12, t_y.to(ureg("N/(m**2)")))
        F_y = q_loc * ds
        Moment_L += F_y * x_loc

    t_y = qs12L[-1] / Wing.ThSpar1 + qs12D[-1] / Wing.ThSpar1
    t_ys12 = np.append(t_ys12, t_y.to(ureg("N/(m**2)")))
    # Moments from section 2 -> 3
    t_xs23 = np.array([])
    t_ys23 = np.array([])
    for i in range(0, len(qs23L)-1):
        q_loc = (qs23L[i] + qs23L[i + 1]) / 2
        s_loc = (s2[i] + s2[i + 1]) / 2
        ds = s2[i + 1] - s2[i]
        x_loc_1, y_loc_1 = Wing.get_xy_from_perim(s2[i]/Wing.Chordlength, Wing.ChSpar1)
        x_loc_2, y_loc_2 = Wing.get_xy_from_perim(s2[i+1]/Wing.Chordlength, Wing.ChSpar1)
        x_loc_1 *= Wing.Chordlength

        x_loc_2 *= Wing.Chordlength

        y_loc_1 *= Wing.Chordlength
        y_loc_2 *= Wing.Chordlength
        Force_angle = np.arctan2(y_loc_2 - y_loc_1, x_loc_2 - x_loc_1)
        x_loc_1 -= x_coor_AC
        x_loc_2 -= x_coor_AC
        F_x = q_loc * ds * np.cos(Force_angle)
        F_y = q_loc * ds * np.sin(Force_angle)
        t_x = (qs23L[i]/Wing.ThSkin) * np.cos(Force_angle) + qs23D[i]/Wing.ThSkin * np.cos(Force_angle)
        t_y = (qs23L[i]/Wing.ThSkin) * np.sin(Force_angle) + qs23D[i]/Wing.ThSkin * np.sin(Force_angle)

        t_xs23 = np.append(t_xs23, t_x.to(ureg("N/(m**2)")))
        t_ys23 = np.append(t_ys23, t_y.to(ureg("N/(m**2)")))

        Moment_L += -F_x*(y_loc_2 + y_loc_1)/2 + F_y*(x_loc_2 + x_loc_1)/2

    t_x = (qs23L[-1] / Wing.ThSkin) * np.cos(Force_angle) + qs23D[-1] / Wing.ThSkin * np.cos(Force_angle)
    t_y = (qs23L[-1] / Wing.ThSkin) * np.sin(Force_angle) + qs23D[-1] / Wing.ThSkin * np.sin(Force_angle)

    t_xs23 = np.append(t_xs23, t_x.to(ureg("N/(m**2)")))
    t_ys23 = np.append(t_ys23, t_y.to(ureg("N/(m**2)")))
    # Moments from section 3->5
    t_ys35 = np.array([])
    for i in range(0, len(qs35L)-1):
        q_loc = (qs35L[i] + qs35L[i + 1]) / 2
        s_loc = (s3[i] + s3[i + 1]) / 2
        ds = s3[i + 1] - s3[i]
        x_loc = Wing.ChSpar2 * Wing.Chordlength - x_coor_AC

        F_y = q_loc * ds

        t_y = qs35L[i] / Wing.ThSpar2 + qs35D[i] / Wing.ThSpar2

        t_ys35 = np.append(t_ys35, t_y.to(ureg("N/(m**2)")))
        Moment_L += -F_y * x_loc
    t_y = qs35L[-1] / Wing.ThSpar2 + qs35D[-1] / Wing.ThSpar2

    t_ys35 = np.append(t_ys35, t_y.to(ureg("N/(m**2)")))
    # Moments from section 5 -> 6
    t_xs56 = np.array([])
    t_ys56 = np.array([])
    for i in range(0, len(qs56L)-1):
        q_loc = (qs56L[i] + qs56L[i + 1]) / 2
        s_loc = (s4[i] + s4[i + 1]) / 2
        ds = s4[i + 1] - s4[i]
        x_loc_1, y_loc_1 = Wing.get_xy_from_perim(s4[i] / Wing.Chordlength, Wing.ChSpar2, reverse=True)
        x_loc_2, y_loc_2 = Wing.get_xy_from_perim(s4[i + 1] / Wing.Chordlength, Wing.ChSpar2, reverse=True)
        x_loc_1 *= Wing.Chordlength
        x_loc_1 -= x_coor_AC
        x_loc_2 *= Wing.Chordlength
        x_loc_2 -= x_coor_AC
        y_loc_1 *= Wing.Chordlength
        y_loc_2 *= Wing.Chordlength

        Force_angle = np.arctan2(y_loc_2 - y_loc_1, x_loc_2 - x_loc_1)
        t_x = qs56L[i] / Wing.ThSkin * np.cos(Force_angle) + qs56D[i] / Wing.ThSkin * np.cos(Force_angle)
        t_y = qs56L[i] / Wing.ThSkin * np.sin(Force_angle) + qs56D[i] / Wing.ThSkin * np.sin(Force_angle)

        t_xs56 = np.append(t_xs56, t_x.to(ureg("N/(m**2)")))
        t_ys56 = np.append(t_ys56, t_y.to(ureg("N/(m**2)")))
        F_x = q_loc * ds * np.cos(Force_angle)
        F_y = q_loc * ds * np.sin(Force_angle)
        Moment_L += -F_x * (y_loc_2 + y_loc_1) / 2 + F_y * (x_loc_2 + x_loc_1) / 2

    t_x = qs56L[-1] / Wing.ThSkin * np.cos(Force_angle) + qs56D[-1] / Wing.ThSkin * np.cos(Force_angle)
    t_y = qs56L[-1] / Wing.ThSkin * np.sin(Force_angle) + qs56D[-1] / Wing.ThSkin * np.sin(Force_angle)

    t_xs56 = np.append(t_xs56, t_x.to(ureg("N/(m**2)")))
    t_ys56 = np.append(t_ys56, t_y.to(ureg("N/(m**2)")))
    # Moment 6 -> 1
    t_ys61 = np.array([])
    for i in range(0, len(qs61L)-1):
        q_loc = (qs61L[i] + qs61L[i + 1]) / 2
        s_loc = (s5[i] + s5[i + 1]) / 2
        ds = s5[i + 1] - s5[i]
        x_loc = Wing.ChSpar1 * Wing.Chordlength - x_coor_AC
        t_y = qs61L[i] / Wing.ThSpar1 + qs61D[i] / Wing.ThSpar1

        t_ys61 = np.append(t_ys61, t_y.to(ureg("N/(m**2)")))
        F_y = q_loc * ds
        Moment_L += F_y * x_loc
    t_y = qs61L[-1] / Wing.ThSpar1 + qs61D[-1] / Wing.ThSpar1
    t_ys61 = np.append(t_ys61, t_y.to(ureg("N/(m**2)")))

    t_xs23 *= ureg("N/(m**2)")
    t_xs56 *= ureg("N/(m**2)")
    t_ys12 *= ureg("N/(m**2)")
    t_ys23 *= ureg("N/(m**2)")
    t_ys35 *= ureg("N/(m**2)")
    t_ys56 *= ureg("N/(m**2)")
    t_ys61 *= ureg("N/(m**2)")
    return Moment_L, t_xs23, t_xs56, t_ys12, t_ys23, t_ys35, t_ys56, t_ys61



Moment_L, t_xs23, t_xs56, t_ys12, t_ys23, t_ys35, t_ys56, t_ys61  = Calc_moment_due_to_shear(s1, s2, s3, s4, s5, qs12L+qs0_L, qs23L+qs0_L, qs35L+qs0_L, qs56L+qs0_L, qs61L+qs0_L, qs12D+qs0_D, qs23D+qs0_D, qs35D+qs0_D, qs56D+qs0_D, qs61D+qs0_D)

print(Moment_L)
# Calculate shear center location   #Tobias
#units checked and correct
def Shear_center(moment_shear):
    shear_center = moment_shear/WingStress.L
    return shear_center

shear_center = Shear_center(Moment_L)
# Calculate Torque                 #Tobias
#units checked and correct
def Torque_for_twist(shear_center):
    T = WingStress.M + WingStress.L * shear_center
    return T

T = Torque_for_twist(shear_center)
# Calculate Rate of Twist          #Tobias
#units checked and correct
def Rate_of_twist(T):
    constant = T/(4*Wing.Area_cell()**2*WingStress.shear_modulus)
    integral = Wing.HSpar1/Wing.ThSpar1
    integral += 2*Wing.length_Skin_x_c(Wing.ChSpar1, Wing.ChSpar2)/Wing.ThSkin
    integral += Wing.HSpar2/Wing.ThSpar2
    dthetadz = constant/integral
    return dthetadz

dthetadz = Rate_of_twist(T)

# Final shear flows in each section
def Final_shaer_flows(qs12L, qs23L, qs35L, qs56L, qs61L, qs12D, qs23D, qs35D, qs56D, qs61D, qs0_L, qs0_D, q_moment):
    qs12 = qs12L + qs12D + qs0_D +qs0_L + q_moment
    qs23 = qs23L + qs23D + qs0_D +qs0_L + q_moment
    qs35 = qs35L + qs35D + qs0_D +qs0_L + q_moment
    qs56 = qs56L + qs56D + qs0_D +qs0_L + q_moment
    qs61 = qs61L + qs61D + qs0_D +qs0_L + q_moment
    return qs12, qs23, qs35, qs56, qs61

qs12, qs23, qs35, qs56, qs61 = Final_shaer_flows(qs12L, qs23L, qs35L, qs56L, qs61L, qs12D, qs23D, qs35D, qs56D, qs61D, qs0_L, qs0_D, q_moment)
# Calculate shear stress

def Get_xy_components(s1, s2, s3, s4, s5, qs12, qs23, qs35, qs56, qs61):

    t_ys12 = np.array([])

    x_coor_AC = 0.25*Wing.Chordlength
    Moment_L = 0
    # Moments from section 1 -> 2
    for i in range(0,len(qs12)-1):
        t_y = qs12[i]/Wing.ThSpar1
        t_ys12 = np.append(t_ys12, t_y.to(ureg("N/(m**2)")))
    t_y = qs12L[-1] / Wing.ThSpar1
    t_ys12 = np.append(t_ys12, t_y.to(ureg("N/(m**2)")))
    # Moments from section 2 -> 3
    t_xs23 = np.array([])
    t_ys23 = np.array([])
    for i in range(0, len(qs23)-1):
        q_loc = (qs23L[i] + qs23L[i + 1]) / 2
        s_loc = (s2[i] + s2[i + 1]) / 2
        ds = s2[i + 1] - s2[i]
        x_loc_1, y_loc_1 = Wing.get_xy_from_perim(s2[i]/Wing.Chordlength, Wing.ChSpar1)
        x_loc_2, y_loc_2 = Wing.get_xy_from_perim(s2[i+1]/Wing.Chordlength, Wing.ChSpar1)
        x_loc_1 *= Wing.Chordlength
        x_loc_2 *= Wing.Chordlength
        y_loc_1 *= Wing.Chordlength
        y_loc_2 *= Wing.Chordlength
        Force_angle = np.arctan2(y_loc_2 - y_loc_1, x_loc_2 - x_loc_1)
        t_x = qs23[i]/Wing.ThSkin * np.cos(Force_angle)
        t_y = qs23[i]/Wing.ThSkin * np.sin(Force_angle)
        t_xs23 = np.append(t_xs23, t_x.to(ureg("N/(m**2)")))
        t_ys23 = np.append(t_ys23, t_y.to(ureg("N/(m**2)")))

    t_x = qs23[-1] / Wing.ThSkin * np.cos(Force_angle)
    t_y = qs23[-1] / Wing.ThSkin * np.sin(Force_angle)

    t_xs23 = np.append(t_xs23, t_x.to(ureg("N/(m**2)")))
    t_ys23 = np.append(t_ys23, t_y.to(ureg("N/(m**2)")))
    # Moments from section 3->5
    t_ys35 = np.array([])
    for i in range(0, len(qs35)-1):
        t_y = qs35[i] / Wing.ThSpar2
        t_ys35 = np.append(t_ys35, t_y.to(ureg("N/(m**2)")))
        Moment_L += -F_y * x_loc
    t_y = qs35L[-1] / Wing.ThSpar2 + qs35D[-1] / Wing.ThSpar2

    t_ys35 = np.append(t_ys35, t_y.to(ureg("N/(m**2)")))
    # Moments from section 5 -> 6
    t_xs56 = np.array([])
    t_ys56 = np.array([])
    for i in range(0, len(qs56)-1):
        q_loc = (qs56[i] + qs56[i + 1]) / 2
        s_loc = (s4[i] + s4[i + 1]) / 2
        ds = s4[i + 1] - s4[i]
        x_loc_1, y_loc_1 = Wing.get_xy_from_perim(s4[i] / Wing.Chordlength, Wing.ChSpar2, reverse=True)
        x_loc_2, y_loc_2 = Wing.get_xy_from_perim(s4[i + 1] / Wing.Chordlength, Wing.ChSpar2, reverse=True)
        x_loc_1 *= Wing.Chordlength
        x_loc_1 -= x_coor_AC
        x_loc_2 *= Wing.Chordlength
        x_loc_2 -= x_coor_AC
        y_loc_1 *= Wing.Chordlength
        y_loc_2 *= Wing.Chordlength

        Force_angle = np.arctan2(y_loc_2 - y_loc_1, x_loc_2 - x_loc_1)
        t_x = qs56[i] / Wing.ThSkin * np.cos(Force_angle)
        t_y = qs56[i] / Wing.ThSkin * np.sin(Force_angle)

        t_xs56 = np.append(t_xs56, t_x.to(ureg("N/(m**2)")))
        t_ys56 = np.append(t_ys56, t_y.to(ureg("N/(m**2)")))

    t_x = qs56L[-1] / Wing.ThSkin * np.cos(Force_angle) + qs56D[-1] / Wing.ThSkin * np.cos(Force_angle)
    t_y = qs56L[-1] / Wing.ThSkin * np.sin(Force_angle) + qs56D[-1] / Wing.ThSkin * np.sin(Force_angle)

    t_xs56 = np.append(t_xs56, t_x.to(ureg("N/(m**2)")))
    t_ys56 = np.append(t_ys56, t_y.to(ureg("N/(m**2)")))
    # Moment 6 -> 1
    t_ys61 = np.array([])
    for i in range(0, len(qs61)-1):
        t_y = qs61[i] / Wing.ThSpar1
        t_ys61 = np.append(t_ys61, t_y.to(ureg("N/(m**2)")))
    t_y = qs61L[-1] / Wing.ThSpar1 + qs61D[-1] / Wing.ThSpar1
    t_ys61 = np.append(t_ys61, t_y.to(ureg("N/(m**2)")))

    t_xs23 *= ureg("N/(m**2)")
    t_xs56 *= ureg("N/(m**2)")
    t_ys12 *= ureg("N/(m**2)")
    t_ys23 *= ureg("N/(m**2)")
    t_ys35 *= ureg("N/(m**2)")
    t_ys56 *= ureg("N/(m**2)")
    t_ys61 *= ureg("N/(m**2)")
    return t_xs23, t_xs56, t_ys12, t_ys23, t_ys35, t_ys56, t_ys61


qs23X, qs56X, qs12Y, qs23Y, qs35Y, qs56Y, qs61Y = Calc_moment_due_to_shear(s1, s2, s3, s4, s5, qs12, qs23, qs35, qs56, qs61)



#Tsia-Wu Failure criterion
def Tsia_Wu(sigma_zs, tau_x, tau_y):
    F11=1/(yield_strength*WingStress.compr_strength)
    F22 = F11
    F12 = -1/2*np.sqrt(F11*F22)
    F1 = 1/(WingStress.yield_strength)-1/(WingStress.compr_strength)
    F2 = 1/(WingStress.yield_strength)-1/(WingStress.compr_strength)
    F44 = 1/WingStress.tau_max**2
    F66 = 1/WingStress.tau_max**2
    sigma1 = sigma_zs
    sigma2 = 0
    sigma3 = 0
    tau12 = tau_y
    tau23 = 0
    tau13 = tau_x
    F = F11 *sigma1**2+F22*(sigma2**2+sigma3**2)+sigma2*sigma3*(2*F22-F44)
    F += 2*F12*sigma1*(sigma3+sigma2)+F1*(sigma1+sigma2) + F2*sigma3
    F += F44*tau23**2 + F66*(tau13**2+tau12**2)
    if F < 1:
        print("No failure occurs")
    else:
        print("Failure occurs")
    return F

F = Tsia_Wu(WingStress.Normal_stress_due_to_bending(0.18, Wing.airfoilordinate(0.18))[0], qs23X[20], qs23Y[20])
print("F =", F)

#plt.plot(s3, qs3)
#plt.show()