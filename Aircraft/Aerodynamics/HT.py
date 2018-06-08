"""         
Name: Horizontal Tail
Department: Aerodynamics
Last updated: 05/06/2018 12:45 by Midas
"""

import sys
import math as m
sys.path.append('../')

toc = 0.09                  # Thickness over chord
x_trl = 0.05                # Transition point lower wing H-tail
x_tru = 0.05                # Transition point upper wing H-tail
C_Nh_alpha = 3.2389         # Normal force coefficient H-tail
Vh_v = 0.925                # V-ratio H-tail
cl_de = 5.34                # Change in C_L due to elevator deflection