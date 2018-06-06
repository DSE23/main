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

def calc_stringer_Inertia(h_str, w_str, t_str):
    
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
    
    I_xy = A_1*
    #I_xy = 
    
    return((I_xx, I_yy))