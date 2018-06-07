# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

Stability derivatives range for level 1 flying qualities
"""

import sys
sys.path.append('../')

from Misc import Init_parm as IP
from Geometry import Wing

# This file will calculate the ranges for l_h,
# Xlemac, S_h and S_v in which StefX will havel level 1 flying qualities

#Input parameters

A = Wing.A
CL_alpha = IP.CL_alpha
MTOW = IP.MTOW
#specific parameters (only check if either Lambda or A changes)


#Stability Derivatives




