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
import scipy as sp

from Geometry import Wing # Import all wing geometry variables

#print(Wing.h_str)

def calc_stringer_inertia(h_str, w_str, t_str):

    # Calculating Area's of stringer components
    A_1 = h_str*t_str
    A_2 = (w_str-t_str)*t_str

    # Calculating Centroid of stringer w.r.t. bottom left corner
    y_bar = (A_1*h_str/2 + A_2*(h_str-0.5*t_str))/(A_1+A_2)
    x_bar = (A_1*0.5*t_str + A_2*(0.5*(w_str - t_str)+t_str))/(A_1+A_2)

    # Calculate I_xx w.r.t. bottom left corner
    I_xx = (t_str * h_str ** 3) / 12 + A_1 * ((h_str / 2)) ** 2 + (
                t_str ** 3 * (w_str - t_str)) / 12 + A_2 * (h_str - 0.5 * t_str) ** 2

    # Calculate I_xx of the centroid of the stringer
    I_xx_centroid = (t_str*h_str**3)/12 + A_1*(((h_str/2)-y_bar))**2 + (t_str**3 * (w_str-t_str))/12 +  A_2*((h_str-0.5*t_str)-y_bar)**2

    I_xx.ito("m**4") # Convert to m^4
    I_xx_centroid.ito("m**4") # Convert to m^4

    # Calculate I_yy of the stringer w.r.t. the bottom left corner
    I_yy = (((w_str - t_str) ** 3) * t_str) / 12 + A_2 * (0.5 * (w_str - t_str) + t_str) ** 2 + (
                t_str ** 3 * h_str) / 12 + A_1 * (0.5 * t_str) ** 2
    # Calculate I_yy of the centroid of the stringer
    I_yy_centroid = (((w_str-t_str)**3) * t_str)/12 + A_2 * (0.5*(w_str - t_str)+t_str - x_bar)**2 + (t_str**3 * h_str)/12 + A_1*(0.5*t_str - x_bar)**2

    I_yy.ito("m**4") # Convert to m^4
    I_yy_centroid.ito("m**4") # Convert to m^4

    # Calculate I_xy of the stringer w.r.t. the bottom-left hand corner
    I_xy = A_1 * (h_str / 2 - y_bar) * (t_str / 2 - 0) + A_2 * ((h_str - 0.5 * t_str) - 0) * (
                (w_str - t_str) * 0.5 + t_str - 0)
    # Calculate I_xy of the stringer w.r.t. the centroid
    I_xy_centroid = A_1 * (h_str / 2 - y_bar) * (t_str / 2 - x_bar) + A_2 * ((h_str - 0.5 * t_str) - y_bar) * ((w_str-t_str)*0.5 + t_str - x_bar)

    I_xy_centroid.ito("m**4") # Convert to m^4

    return ((I_xx, I_yy, I_xy), (I_xx_centroid, I_yy_centroid, I_xy_centroid))

def Calc_spar_inertia(HSpar,TSpar,Centroid,ChSpar,zs):           #Input height spar (m), thickness spar (m), location centroid w.r.t. chord and chordwise location (-) spar respectively (-)
    Ixx = (1/12)*TSpar*(HSpar**3)             #Calculation of Ixx
    Iyy = (1/12)*HSpar*(TSpar**3)               #Caclulation of Iyy w/o steiner term
    Iyysteiner = TSpar*HSpar*(abs((Centroid*length_chord(zs))-(ChSpar*length_chord(zs)))**2)      #Calculation of steiner term Iyy
    Iyy = Iyy + Iyysteiner                              #Adding both Iyy moments of inertia together
    return Ixx, Iyy

def Calc_skin_inertia_Ixx(Spar1, Spar2):
    n = 100 #number of sections
    dx = ((Spar2-Spar1)/n)
    x = Spar1
    Ixx = 0
    for i in range(n):
        x = x + dx
        dxlength=dx*Wing.Chordlength
        y = ((Wing.airfoilordinate(x-dx)+Wing.airfoilordinate(x))/2)*Wing.Chordlength
        dy = abs(Wing.airfoilordinate(x-dx)-Wing.airfoilordinate(x))*Wing.Chordlength
        length = np.sqrt(dxlength**2+dy**2)
        dIxx = length*Wing.ThSkin*(y**2)
        Ixx = Ixx+dIxx
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
        dIyy = length * Wing.ThSkin*((abs(xlength-(Wing.centroid*Wing.Chordlength)))**2)
        Iyy = Iyy + dIyy

    return Iyy




# returns stiffener x,y locations and rotation
# return z_y_angle_coords  # [(stringer0 z,y,rot),(stringer1 x,y,rot)] m,m,rad
def stif_loc(z, t_sk, n_st, x):
    total_perimeter = sp.integrate.quad(Wing.airfoilordinate(x), Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), Wing.Chord_loc_Spar(z,Wing.Spar2R,Wing.Spar2T)) #m

    spacing = total_perimeter / ((n_st + 1) / 2)
    x_y_angle_coords = []
    for i in range(n_st):
        local_spacing = i * spacing
        rot_angle = Wing.Angle(x) + radians(180)
        x_coordinate = Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T) + sp.integrate.quad(math.cos(Wing.angle(x)),Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), localspacing) 
        #x_coordinate = (-1) * (local_spacing - circle_perim) * cos(atan(0.5 * h / (C_a - 0.5 * h)))
        y_coordinate = Wing.airfoilordinate(Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T)+local_spacing)-sp.integrate.quad(math.sin(Wing.angle(x)),Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), localspacing)
        
        apnd_itm = (x_coordinate, y_coordinate, rot_angle)
        x_y_angle_coords.append(apnd_itm)
        apnd_itm = (x_coordinate, -y_coordinate, -rot_angle)
        x_y_angle_coords.append(apnd_itm)

        # print "Stif.", i, "\t x:", x_coordinate, "\t y:", y_coordinate, "\t angle:", degrees(rot_angle)

    return x_y_angle_coords  # [(stringer0 x,y,rot),(stringer1 x,y,rot), ...]

print('this is', Calc_skin_inertia_Ixx(Wing.ChSpar1,Wing.ChSpar2))
print('this is', Calc_skin_inertia_Iyy(Wing.ChSpar1,Wing.ChSpar2))
#print(calc_stringer_Inertia(Q_("50 mm"),Q_("20 mm"),Q_("2 mm")))

# VERIFICATION TESTS
class InertiaTestCase(unittest.TestCase):
    def setUp(self):
        self.Inertias = calc_stringer_inertia(Q_("50 mm"),Q_("20 mm"),Q_("2 mm"))
        self.I_xx_centr = self.Inertias[1][0]
        self.I_xx_centr.ito(ureg("mm**4"))
        self.I_yy_centr = self.Inertias[1][1]
        self.I_yy_centr.ito(ureg("mm**4"))
        self.I_xy_centr = self.Inertias[1][2]
        self.I_xy_centr.ito(ureg("mm**4"))

    def test_Ixx_centroid(self):
        self.assertAlmostEqual(self.I_xx_centr.magnitude, 36092.39, 1,
                         msg="Verification of I_xx for stringer with SolidWorks FAILED")
    def test_Iyy_centroid(self):
        self.assertAlmostEqual(self.I_yy_centr.magnitude, 3652.39, 1,
                               msg="Verification of I_yy for stringer with SolidWorks FAILED")
    def test_Ixy_centroid(self):
        self.assertAlmostEqual(self.I_xy_centr.magnitude, 6352.94, 1,
                               msg="Verification of I_xy for stringer with SolidWorks FAILED")

if __name__ == '__main__':
    unittest.main()