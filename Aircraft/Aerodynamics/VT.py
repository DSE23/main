"""                  
Name: Vertical Tail
Department: Aerodynamics
Last updated: 05/06/2018 12:45 by Midas
"""

import sys
import math as m
sys.path.append('../')

toc = 0.15                  # Thickness over chord
x_trl = 0.05                # Transition point lower wing V-tail
x_tru = 0.05                # Transition point upper wing V-tail
C_Yv_alpha = 2.106          # Side force coefficient V-tail
Vv_V = 0.925                # V-ratio V-tail
q_dsigma_dbeta = 0.8736     # Sidewash gradient (1-dsigma/dbeta)*q/q_infty
cy_dr = 5.5714              # Change in Cy due to rudder deflection
