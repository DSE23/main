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

# Calculate Boom Area's             #Midas

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