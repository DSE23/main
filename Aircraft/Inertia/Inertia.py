# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

Inertia parameters

"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import math as m
import numpy as np

I_yy = Q_("1594 kg m**2")
I_xx = Q_("1089 kg m**2")
I_zz = Q_("2629 kg m**2")
