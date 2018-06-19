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
# Calculate Boom Area's

# Calculate base shear flow for every section of the wing box

# Calculate correcting shear flow (qs0)

# Add correcting shear flow to base shear flows

# Compute moments around a.c. caused by shear forces due to shear flows

# Calculate shear center location
def Shear_center(moment_shear,):
    shear_center = moment_shear/WingStress.L
    return shear_center

# Calculate Torque
def Torque_for_twist(shear_center):
    T = WingStress.M + WingStress.L * shear_center
    return T

# Calculate Rate of Twist
def Rate_of_twist(T):
    constant = T/(4*Wing.Area_cell()*WingStress.shear_modulus)
    integral = Wing.HSpar1/Wing.ThSpar1
    integral += 2*Wing.length_Skin_x_c/Wing.ThSkin
    integral += Wing.HSpar2/Wing.ThSpar2
    dthetadz = constant/integral
    return dthetadz
    
# Calculate shear stress