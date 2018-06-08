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
CL_max_hld = CL_max + dCL_flaps + dCL_slats
alpha_stall = 18.20                     # Stall angle of attack
CL_alpha = 4.6532                       # CL change due to AoA change
Oswald_e = 0.792                        # Oswald efficiency factor
C_Nw_alpha = 4.6532                     # Normal force coef with respect to alpha
de_da = 0.806                           # Downwash gradient
cl_da = 3.9282                          # Change in CL due to Aileron deflection

