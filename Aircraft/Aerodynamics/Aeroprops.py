# -*- coding: utf-8 -*-
"""
Name: Aerodynamic Props
Department: Aero
Last updated: 11/06/2018 11.21 by Emma
"""


import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import math as m
import numpy as np


q_qinf_ratio = 0.8558 #ratio dynamic pressure ht over infinity dyn pressure
CL_alpha_wing = Q_("4.65 1/rad") #slope wing lift co
CL_alpha_ht = Q_("3.23 1 / rad") #slope ht lift co, downwash not included
de_da = 0.806 #downwash angle change over angle of attack change
CD0_wing = Q_("0.0075473919026861095")
CD0_ht = Q_("0.001798445721745142")
CD0_vt = Q_("0.0007209779478792948")
CD_canopy = Q_("0.04 dimensionless")
CD_lg = Q_("0.032572960000000005 dimensionless")
CD0_tot = Q_("0.015547346906243672 dimensionless")






























