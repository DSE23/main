"""                  
Name: Wing
Department: Aerodynamics
Last updated: 05/06/2018 12:45 by Midas
"""

import sys
import math as m
sys.path.append('../')

CL_max = 1.3454                         # CL max clean config
dCL_flaps = 0.567                       # Delta CL due to flaps
dCL_slats = 0.64                        # Delta CL due to the slats
alpha_stall = 18.20                     # Stall angle of attack
CL_alpha = 4.6532                       # CL change due to AoA change
CD_0 = 0.028                            # Zero lift drag coeff
