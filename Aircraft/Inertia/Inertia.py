"""                  
Name: Inertia
Department: Inertia
Last updated: 06/06/2018 12:41 by Midas
"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import numpy as np
import scipy.integrate as integrate

from Geometry import Wing # Import all wing geometry variables

#print(Wing.h_str)


def calc_stringer_inertia(h_str, w_str, t_str):

    # Calculating Area's of stringer components
    A_1 = h_str*t_str
    A_2 = (w_str-t_str)*t_str
    
    # Calculating Centroid of stringer w.r.t. bottom left corner
    y_bar = (A_1*h_str/2 + A_2*(h_str-0.5*t_str))/(A_1+A_2)
    x_bar = (A_1*0.5*t_str + A_2*(w_str/2 + t_str))/(A_1+A_2)

    
    I_xx = (t_str*h_str**3)/3 # Ixx inertia of vertical part of the stringer, referenced to the bottom
    I_xx += t_str*w_str*h_str**2 # Ixx inertia of horizontal part with thin walled approximation
    
    I_yy = (w_str**3*t_str)/3 # Iyy Inertia of horizontal part referenced to bottom-left corner
    I_yy += 0 # Iyy of vertical part is negligible (thin-walled approximation)
    
    I_xy = A_1

    print(y_bar)
    print(x_bar)
    
    I_xx = (t_str*h_str**3)/12 # Ixx inertia of vertical part of the stringer
    I_xx += A_1*h_str**2
    I_xx += (t_str*(w_str-t_str))*(h_str-0.5*t_str)**2 # Ixx inertia of horizontal part with thin walled approximation
    I_xx.ito("m**4")
    
    I_yy = (((w_str-t_str)**3) * t_str)/12 + A_2 * (0.5*(w_str))**2
    # I_yy = ((w_str-t_str)**3 * t_str)/3 # Iyy Inertia of horizontal part referenced to bottom-left corner
    # I_yy += A_1*(0.5*t_str)**2 # Iyy of vertical part is negligible (thin-walled approximation)
    I_yy.ito("m**4")
    
    I_xy = A_1*(h_str/2 - 0)*(t_str/2 - 0)
    I_xy += A_2*((h_str-0.5*t_str) - 0)*((t_str+0.5*w_str) -0)
    I_xy.ito("m**4")
 
 
    
    return((I_xx, I_yy, I_xy))
    

    return I_xx, I_yy, I_xy


 
    


def Calc_skin_inertia(Spar1, Spar2):
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

print('this is', Calc_skin_inertia(Wing.ChSpar1,Wing.ChSpar2))
#print(calc_stringer_Inertia(Q_("50 mm"),Q_("20 mm"),Q_("2 mm")))




# returns stiffener x,y locations and rotation
# return z_y_angle_coords  # [(stringer0 z,y,rot),(stringer1 x,y,rot)] m,m,rad
def stif_loc(z, t_sk, n_st):
    total_perimeter = integrate.quad(Wing.airfoilordinate(x), Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), Wing.Chord_loc_Spar(z,Wing.Spar2R,Wing.Spar2T)) #m

    spacing = total_perimeter / ((n_st + 1) / 2)
    x_y_angle_coords = []
    for i in xrange(6):
        local_spacing = i * spacing
        if local_spacing < circle_perim:
            angle = (local_spacing / circle_perim) * radians(90)
            x_coordinate = -1 * (0.5 * h - (0.5 * h - t_sk + cos(angle) * (0.5 * h - t_sk)))
            y_coordinate = sin(angle) * (0.5 * h - t_sk)
            rot_angle = angle + radians(90)

        else:
            rot_angle = atan(0.5 * h / (C_a - 0.5 * h)) - radians(180)
            x_coordinate = (-1) * (local_spacing - circle_perim) * cos(atan(0.5 * h / (C_a - 0.5 * h)))
            y_coordinate = h / 2 - (local_spacing - circle_perim) * sin(atan(0.5 * h / (C_a - 0.5 * h)))

        apnd_itm = (x_coordinate, y_coordinate, rot_angle)
        x_y_angle_coords.append(apnd_itm)
        if i > 0:
            apnd_itm = (x_coordinate, -y_coordinate, -rot_angle)
            x_y_angle_coords.append(apnd_itm)

        # print "Stif.", i, "\t x:", x_coordinate, "\t y:", y_coordinate, "\t angle:", degrees(rot_angle)

    return x_y_angle_coords  # [(stringer0 x,y,rot),(stringer1 x,y,rot), ...]


