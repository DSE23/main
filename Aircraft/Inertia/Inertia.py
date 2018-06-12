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
from Geometry import Geometry

I_yy = Q_("1594 kg m**2")
I_xx = Q_("1089 kg m**2")
I_zz = Q_("2629 kg m**2")

#New inertia Calculation

#Fuselage

Afusi = ([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
y_fus1 = -1
y_fus2 = 0
z_fus1 = -1
z_fus2 = 0
z_fus3 = 1
SumAfus = sum(Afusi)
W_fus = Geometry.Masses.W_fus
W_ifus = W_fus * Afusi/SumAfus
m_ifus = W_ifus/8

theta1 = np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))
theta2 = m.radians(90) - theta1
theta3 = np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))
theta4 = m.radians(90) - theta3
theta5 = theta4
theta6 = theta3
theta7 = theta2
theta8 = theta1

alpha1 = theta1
alpha2 = (theta1+theta2)/2
alpha3 = alpha2
alpha4 = (theta3 + theta4)/2
alpha5 = theta4
alpha6 = alpha4
alpha7 = alpha3
alpha8 = alpha2

