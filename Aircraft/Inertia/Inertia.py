"""                  
Name: Inertia
Department: Inertia
Last updated: 05/06/2018 12:41 by Midas
"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder

from Geometry import Wing # Import all wing geometry variables

print(Wing.h_str)


def calc_stringer_inertia(h_str, w_str, t_str):

    # Calculating Area's of stringer components
    A_1 = h_str*t_str
    A_2 = (w_str-t_str)*t_str
    print(A_2)
    # Calculating Centroid of stringer w.r.t. bottom left corner
    y_bar = (A_1*h_str/2 + A_2*(h_str-0.5*t_str))/(A_1+A_2)
    x_bar = (A_1*0.5*t_str + A_2*(w_str/2 + t_str))/(A_1+A_2)
    
    
    
    I_xx = (t_str*h_str**3)/12 # Ixx inertia of vertical part of the stringer
    I_xx += A_1*(h_str/2)**2
    I_xx += (t_str**3 * (w_str-t_str))/12  # Ixx inertia of horizontal part with thin walled approximation
    I_xx += A_2*(h_str+0.5*t_str)**2
    
    I_xx.ito("m**4")
    
    I_yy = (((w_str-t_str)**3) * t_str)/12 + A_2 * (0.5*(w_str))**2
    # I_yy = ((w_str-t_str)**3 * t_str)/3 # Iyy Inertia of horizontal part referenced to bottom-left corner
    # I_yy += A_1*(0.5*t_str)**2 # Iyy of vertical part is negligible (thin-walled approximation)


    I_yy.ito("m**4")

    I_xy = A_1 * (h_str / 2 - 0) * (t_str / 2 - 0)
    I_xy += A_2 * ((h_str - 0.5 * t_str) - 0) * ((t_str + 0.5 * w_str) - 0)
    I_xy.ito("m**4")

    return ((I_xx, I_yy, I_xy))

def Calc_spar_inertia(HSpar,TSpar,Centroid,ChSpar,zs):           #Input height spar (m), thickness spar (m), location centroid w.r.t. chord and chordwise location (-) spar respectively (-)
    Ixx = (1/12)*TSpar*(HSpar**3)             #Calculation of Ixx
    Iyy = (1/12)*HSpar*(TSpar**3)               #Caclulation of Iyy w/o steiner term
    Iyysteiner = TSpar*HSpar*(abs((Centroid*length_chord(zs))-(ChSpar*length_chord(zs)))**2)      #Calculation of steiner term Iyy
    Iyy = Iyy + Iyysteiner                              #Adding both Iyy moments of inertia together
    return Ixx, Iyy

print(calc_stringer_Inertia(Q_("50 mm"),Q_("20 mm"),Q_("2 mm")))
