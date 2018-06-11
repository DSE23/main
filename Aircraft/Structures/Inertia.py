
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

def get_xy_from_perim(perim_val, start_x=0, dat_file_name="../Airfoil.dat"):
    """
        NOTE: Special function for Tobias
        This function returns the x and y coordinates for a given perimeter


        Arguments: -perim_val (x coordinate)
                   -start_x: x_value from where the perimeter integration needs to start
                   -dat_file_name (name of the airfoil data file)

        Returns: perim (Perimiter value)
    """

    Air_data = np.genfromtxt(dat_file_name)  # Imports datapoints from airfoil data file

    x_coords = Air_data[:81, 0]  # We only care about 1 half of the airfoil
    x_coords = np.flip(x_coords, 0)  # Flip them so they are in a good order
    y_coords = Air_data[:81, 1]  # We only care about 1 half of the airfoil
    y_coords = np.flip(y_coords, 0)  # Flip them so they are in a good order
    p = interp1d(x_coords, y_coords, kind='cubic')  # Generate a poly spline based on the airfoil points

    perim = 0  # Set initial perimiter size to 0
    step = 0.0001  # Step size for algorithm: increase will lead to faster computing times
    min_val = 10
    for x_c in np.arange(start_x, 1, step):
        perim += np.sqrt((step) ** 2 + (p(x_c + step) - p(x_c)) ** 2)
        if abs(perim - perim_val) < min_val:
            x_coord = 0.5 * (x_c + x_c + step)
            y_coord = 0.5 * (p(x_c) + p(x_c + step))
            min_val = abs(perim - perim_val)
        else:
            break

    return(x_coord, y_coord)

def get_perim_from_x(x_coor, dat_file_name="../Airfoil.dat"):
    """
    This function returns the perimeter value from the LE until the specified x-coordinate

    Arguments: x_coor (x coordinate)
               dat_file_name (name of the airfoil data file)

    Returns: perim (Perimiter value)
    """
    Air_data = np.genfromtxt(dat_file_name)  # Imports datapoints from airfoil data file

    x_coords = Air_data[:81, 0]  # We only care about 1 half of the airfoil
    x_coords = np.flip(x_coords, 0)  # Flip them so they are in a good order
    y_coords = Air_data[:81, 1]  # We only care about 1 half of the airfoil
    y_coords = np.flip(y_coords, 0)  # Flip them so they are in a good order
    p = interp1d(x_coords, y_coords, kind='cubic')  # Generate a poly spline based on the airfoil points

    perim = 0  # Set initial perimiter size to 0
    step = 0.0001  # Step size for algorithm: increase will lead to faster computing times
    # but lower accuracy
    # This for loop calculates the perimiter until the specified x-coordinate
    for x_c in np.arange(0.0, x_coor, step):
        perim += np.sqrt((step) ** 2 + (p(x_c + step) - p(x_c)) ** 2)

    return perim


def get_coord_from_perim(n_st, start_x, end_x, chord_l, dat_file_name="../Airfoil.dat"):
    """
    This function returns list of coordinate values where a stiffener is placed
    based on the spar locations and number of stiffeners. The stiffeners will
    be distributed equally along the perimiter section between the spars.

    Arguments: -n_st: number of stiffeners
               -start_x: x location of the LE spar
               -end_X: x location of the TE spar
               -dat_file_name: name of the airfoil data file

    Returns: List of x and y coordinates where the stiffeners are placed
    """
    Air_data = np.genfromtxt(dat_file_name)


    x_coords = Air_data[:81, 0]
    x_coords = np.flip(x_coords, 0)
    y_coords = Air_data[:81, 1]
    y_coords = np.flip(y_coords, 0)

    final_x_coords = np.array([])
    final_y_coords = np.array([])
    final_angles = np.array([])
    p = interp1d(x_coords, y_coords, kind='cubic')

    final_x_y_angle_coord = np.array([])
    perim = 0
    stif_perim = get_perim_from_x(end_x) - get_perim_from_x(start_x)
    perim_spacing = stif_perim / (n_st + 1)
    step = 0.0001
    i = 0
    for x_c in np.arange(start_x, end_x, step):
        perim += np.sqrt((step) ** 2 + (p(x_c + step) - p(x_c)) ** 2)
        if i >= n_st:
            break
        if abs(perim - perim_spacing) < 10e-5:
            x_coord = 0.5 * (x_c + x_c + step)
            y_coord = 0.5 * (p(x_c) + p(x_c + step))
            slope_angle = -np.arctan((p(x_c + step) - p(x_c)) / (step)) + np.pi
            final_x_coords = np.append(final_x_coords, x_coord)
            final_y_coords = np.append(final_y_coords, y_coord)
            final_angles = np.append(final_angles, slope_angle)
            perim = 0
            i += 1

    # print(perim)
    a_ran = np.arange(0, 1, 0.0001)
    plt.plot([start_x, start_x], [0, p(start_x)], 'b')
    plt.plot([end_x, end_x], [0, p(end_x)], 'b')
    plt.plot(x_coords, y_coords, '+')
    plt.plot(a_ran, p(a_ran), 'r')
    plt.plot(final_x_coords, final_y_coords, 'go')
    plt.axis((0, 1, 0, 1))

    plt.show()
    return (final_x_coords*chord_l, final_y_coords*chord_l, final_angles)

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
    I_xy_centroid = A_1 * (-(0.5*(b_1) + 0.5*t_str))*(0.5*h_1 - y_bar) + A_3*(0.5*(b_1) + 0.5*t_str)*((h_str-0.5*h_3) - y_bar)

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

#def stiffeners_centroid(x_y_angle_coords, h_st, w_st, t_st):




def calc_total_stringer_inertia(x_y_angle_coords, stringer_inertias):

    x_coords = x_y_angle_coords[0]
    y_coords = x_y_angle_coords[1]
    angles = x_y_angle_coords[2]

    I_xx_stringer = stringer_inertias[0][0]
    I_yy_stringer = stringer_inertias[0][1]
    I_xy_stringer = stringer_inertias[0][2]
    stringer_area = stringer_inertias[2]
    I_XX_TOT = Q_("0 m**4")
    I_YY_TOT = Q_("0 m**4")
    for i in range(len(x_coords)):
        loc_I_xx, loc_I_yy, loc_I_xy = axis_transformation(I_xx_stringer, I_yy_stringer, I_xy_stringer, angles[i])
        I_XX_TOT += loc_I_xx + stringer_area*(y_coords[i])**2
        #I_YY_TOT += loc_I_yy + stringer_area*()

    return I_XX_TOT
# returns stiffener x,y locations and rotation
# return z_y_angle_coords  # [(stringer0 z,y,rot),(stringer1 x,y,rot)] m,m,rad

print(calc_total_stringer_inertia(get_coord_from_perim(5, 0.2, 0.6, Q_("7 m")), calc_stringer_inertia(Q_("30 mm"), Q_("30 mm"), Q_("2 mm"))))

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

    def test_Ixx_centroid(self):
        self.assertAlmostEqual(self.I_xx_centr.magnitude, 15494.67, 1,
                               msg="Verification of I_xx for stringer with SolidWorks FAILED")

    def test_Iyy_centroid(self):
        self.assertAlmostEqual(self.I_yy_centr.magnitude, 4518.67, 1,
                               msg="Verification of I_yy for stringer with SolidWorks FAILED")

    def test_Ixy_centroid(self):
        self.assertAlmostEqual(self.I_xy_centr.magnitude,  6272.00, 1,
                               msg="Verification of I_xy for stringer with SolidWorks FAILED")


if __name__ == '__main__':
    unittest.main()