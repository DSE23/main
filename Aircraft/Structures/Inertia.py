
"""
Name: Inertia
Department: Inertia
Last updated: 06/06/2018 12:41 by Midas
"""
import sys
import unittest
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import numpy as np
from scipy import integrate
from Structures import Wing
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

# print(Wing.h_str)
Iyy_aircraft = Q_("1492.8 kg/m/m")
Ixx_aircraft = Q_("1016.9 kg/m/m")
Izz_aircraft = Q_("2447.2 kg/m/m")

def calc_stringer_inertia(h_str, w_str, t_str):

    b_1 = w_str/2 - 0.5*t_str
    h_1 = t_str

    b_2 = t_str
    h_2 = h_str

    b_3 = b_1
    h_3 = h_1

    # Calculating Area's of stringer components
    A_1 = b_1*h_1
    A_2 = b_2*h_2
    A_3 = A_1

    # Calculating Centroid of stringer w.r.t. bottom corner
    y_bar = h_str/2
    x_bar = 0

    # Calculate I_xx w.r.t. bottom corner
    I_xx = (b_1*h_1**3)/12 + A_1*(h_1/2)**2 + (b_2*h_2**3)/12 + A_2*(h_2/2)**2 + (b_3*h_3**3)/12 + A_3*(h_str-0.5*h_3)**2

    # Calculate I_xx of the centroid of the stringer
    I_xx_centroid = (b_1*h_1**3)/12 + A_1*(h_1/2 - y_bar)**2 + (b_2*h_2**3)/12 + A_2*(h_2/2 - y_bar)**2 + (b_3*h_3**3)/12 + A_3*(h_str-0.5*h_3 - y_bar)**2

    I_xx.ito("m**4")  # Convert to m^4
    I_xx_centroid.ito("m**4")  # Convert to m^4

    # Calculate I_yy of the stringer w.r.t. the bottom left corner
    I_yy = (b_1**3*h_1)/12 + A_1*(b_1*0.5 + 0.5*t_str)**2 + (b_2**3*h_2)/12 + (b_3**3*h_3)/12 + A_3*(b_3*0.5+0.5*t_str)**2
    # Calculate I_yy of the centroid of the stringer
    I_yy_centroid = (b_1**3*h_1)/12 + A_1*(b_1*0.5 + 0.5*t_str)**2 + (b_2**3*h_2)/12 + (b_3**3*h_3)/12 + A_3*(b_3*0.5+0.5*t_str)**2

    I_yy.ito("m**4")  # Convert to m^4
    I_yy_centroid.ito("m**4")  # Convert to m^4

    # Calculate I_xy of the stringer w.r.t. the bottom-left hand corner
    I_xy = A_1 * (-(0.5*(b_1) + 0.5*t_str))*(0.5*h_1 - 0) + A_3*(0.5*(b_1) + 0.5*t_str)*((h_str-0.5*h_3) - 0)
    # Calculate I_xy of the stringer w.r.t. the centroid
    I_xy_centroid = A_1 * (-(0.5*(b_1) + 0.5*t_str) - 0)*(0.5*h_1 - y_bar) + A_3*((0.5*(b_1) + 0.5*t_str) - 0)*((h_str-0.5*h_3) - y_bar)

    I_xy_centroid.ito("m**4")  # Convert to m^4

    return ((I_xx, I_yy, I_xy), (I_xx_centroid, I_yy_centroid, I_xy_centroid), A_1+A_2+A_3)


def Calc_spar_inertia(HSpar, TSpar, ChSpar,
                      zs):  # Input height spar (m), thickness spar (m), location centroid w.r.t. chord and chordwise location (-) spar respectively (-)
    Ixx = (1 / 12) * TSpar * (HSpar ** 3)  # Calculation of Ixx
    Iyy = (1 / 12) * HSpar * (TSpar ** 3)  # Caclulation of Iyy w/o steiner term
    Iyysteiner = TSpar * HSpar * (abs((Wing.centroid * Wing.length_chord(zs)) - (
    ChSpar * Wing.length_chord(zs))) ** 2)  # Calculation of steiner term Iyy
    Iyy = Iyy + Iyysteiner  # Adding both Iyy moments of inertia together
    return Ixx, Iyy


def Calc_skin_inertia_Ixx(Spar1, Spar2):
    n = 100  # number of sections
    dx = ((Spar2 - Spar1) / n)
    x = Spar1
    Ixx = 0
    for i in range(n):
        x = x + dx
        dxlength = dx * Wing.Chordlength
        y = ((Wing.airfoilordinate(x - dx) + Wing.airfoilordinate(x)) / 2) * Wing.Chordlength
        dy = abs(Wing.airfoilordinate(x - dx) - Wing.airfoilordinate(x)) * Wing.Chordlength
        length = np.sqrt(dxlength ** 2 + dy ** 2)
        dIxx = length * Wing.ThSkin * (y ** 2)
        Ixx = Ixx + dIxx
    Ixx = Ixx * 2
    return Ixx


def Calc_skin_inertia_Iyy(Spar1, Spar2):
    n = 100  # number of sections
    dx = ((Spar2 - Spar1) / n)
    x = Spar1
    Iyy = 0
    for i in range(n):
        x = x + dx
        xlength = x * Wing.Chordlength
        dxlength = dx * Wing.Chordlength
        y = ((Wing.airfoilordinate(x - dx) + Wing.airfoilordinate(x)) / 2) * Wing.Chordlength
        dy = abs(Wing.airfoilordinate(x - dx) - Wing.airfoilordinate(x)) * Wing.Chordlength
        length = np.sqrt(dxlength ** 2 + dy ** 2)
        dIyy = length * Wing.ThSkin * ((abs(xlength - (Wing.centroid * Wing.Chordlength))) ** 2)
        Iyy = Iyy + dIyy
    Iyy = Iyy * 2
    return Iyy

# function for transforming axis to a rotated version
# input MMOI I_zz, I_yy, I_zy, and rotation angle rot_angle
# outputs new rotated I_uu, I_vv, I_uv
def axis_transformation(I_xx, I_yy, I_xy, rot_angle):
    # Axis transformation for rotated axis system used for Inertia calculations
    I_uu = (I_xx + I_yy) * 0.5 + (I_xx - I_yy) * 0.5 * np.cos(2 * rot_angle) - I_xy * np.sin(2 * rot_angle)
    I_vv = (I_xx + I_yy) * 0.5 - (I_xx - I_yy) * 0.5 * np.cos(2 * rot_angle) + I_xy * np.sin(2 * rot_angle)
    I_uv = (I_xx - I_yy) * 0.5 * np.sin(2 * rot_angle) + I_xy * np.cos(2 * rot_angle)
    return I_uu, I_vv, I_uv

def calc_total_stringer_inertia(x_y_angle_coords, stringer_inertias, h_str, w_str, t_str):

    x_coords = x_y_angle_coords[0]
    x_coords = np.append(x_coords, x_coords)
    x_coords *= ureg.meter
    print(x_coords)
    y_coords = x_y_angle_coords[1]
    y_coords = np.append(y_coords, -y_coords)
    y_coords *= ureg.meter
    angles = x_y_angle_coords[2]
    angles = np.append(angles, angles)
    angles *= ureg("rad")

    I_xx_stringer = stringer_inertias[0][0]
    I_yy_stringer = stringer_inertias[0][1]
    I_xy_stringer = stringer_inertias[0][2]
    stringer_area = stringer_inertias[2]
    I_XX_TOT = Q_("0 m**4")
    I_YY_TOT = Q_("0 m**4")
    I_XY_TOT = Q_("0 m**4")
    Y_CEN = Q_("0 m")
    for i in range(len(x_coords)):
        loc_I_xx, loc_I_yy, loc_I_xy = axis_transformation(I_xx_stringer, I_yy_stringer, I_xy_stringer, angles[i])
        I_XX_TOT += loc_I_xx + stringer_area*(y_coords[i])**2
        I_YY_TOT += loc_I_yy + stringer_area*(x_coords[i] - Wing.centroid*Wing.Chordlength)**2
        I_XY_TOT += loc_I_xy + stringer_area*(x_coords[i] - Wing.centroid*Wing.Chordlength)*(y_coords[i] - Y_CEN)

    return I_XX_TOT, I_YY_TOT, I_XY_TOT
# returns stiffener x,y locations and rotation
# return z_y_angle_coords  # [(stringer0 z,y,rot),(stringer1 x,y,rot)] m,m,rad

I_XX_TOT_str, I_YY_TOT_str, I_XY_TOT_str = calc_total_stringer_inertia(Wing.get_coord_from_perim(5, 0.2, 0.8, Q_("7 m")), calc_stringer_inertia(Q_("30 mm"), Q_("30 mm"), Q_("2 mm")), Q_("30 mm"), Q_("30 mm"), Q_("2 mm"))
I_XX_Spar1, I_YY_Spar1 = Calc_spar_inertia(Wing.HSpar1, Wing.ThSpar1, Wing.ChSpar1, Wing.z)
I_XX_Spar2, I_YY_Spar2 = Calc_spar_inertia(Wing.HSpar2, Wing.ThSpar2, Wing.ChSpar2, Wing.z)

I_XX_Skin = Calc_skin_inertia_Ixx(Wing.ChSpar1, Wing.ChSpar2)
I_YY_Skin = Calc_skin_inertia_Iyy(Wing.ChSpar1, Wing.ChSpar2)

I_XX_TOTAL = I_XX_TOT_str + I_XX_Spar1 + I_XX_Spar2 + I_XX_Skin
I_YY_TOTAL = I_YY_TOT_str + I_YY_Spar1 + I_YY_Spar2 + I_YY_Skin

print("I_XX TOTAL:", I_XX_TOTAL)
print("I_YY TOTAL:", I_YY_TOTAL)
#print('this is', Calc_skin_inertia_Ixx(Wing.ChSpar1, Wing.ChSpar2))
#print('this is', Calc_skin_inertia_Iyy(Wing.ChSpar1, Wing.ChSpar2))


# print(calc_stringer_Inertia(Q_("50 mm"),Q_("20 mm"),Q_("2 mm")))
#stiff_x_y_coord = get_coord_from_perim(5, 0.2, 0.6, Q_("7 m"))
#print(stiff_x_y_coord)


#VERIFICATION TESTS
class InertiaTestCase(unittest.TestCase):
    def setUp(self):
        self.Inertias = calc_stringer_inertia(Q_("30 mm"), Q_("30 mm"), Q_("2 mm"))
        self.I_xx_centr = self.Inertias[1][0]
        self.I_xx_centr.ito(ureg("mm**4"))
        self.I_yy_centr = self.Inertias[1][1]
        self.I_yy_centr.ito(ureg("mm**4"))
        self.I_xy_centr = self.Inertias[1][2]
        self.I_xy_centr.ito(ureg("mm**4"))

        self.I_xx_transf, self.I_yy_transf, self.I_xy_transf = axis_transformation(15494.67, 4518.67,6272 , -30*np.pi/180)


    def test_Ixx_centroid(self):
        self.assertAlmostEqual(self.I_xx_centr.magnitude, 15494.67, 1,
                               msg="Verification of I_xx for stringer with SolidWorks FAILED")

    def test_Iyy_centroid(self):
        self.assertAlmostEqual(self.I_yy_centr.magnitude, 4518.67, 1,
                               msg="Verification of I_yy for stringer with SolidWorks FAILED")

    def test_Ixy_centroid(self):
        self.assertAlmostEqual(self.I_xy_centr.magnitude,  6272.00, 1,
                               msg="Verification of I_xy for stringer with SolidWorks FAILED")

    def test_Transformation(self):
        self.assertAlmostEqual(self.I_xx_transf, 18182.38, 1,
                               msg="Verification of transformation formula with SolidWorks FAILED")
        self.assertAlmostEqual(self.I_yy_transf, 1830.96, 1,
                               msg="Verification of transformation formula with SolidWorks FAILED")
        self.assertAlmostEqual(self.I_xy_transf, -1616.75, 1,
                               msg="Verification of transformation formula with SolidWorks FAILED")



if __name__ == '__main__':
    unittest.main()