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

# New inertia Calculation

# Fuselage

x_fus = np.linspace(0, 6.2, 16)
y_fus1 = [-400, -502, -520, -520, -509, -482, -439, -390, -342, -292,
          -243, -194, -145, -96, -47, -7.3] * Q_("mm")  # left fuselage y-coord
y_fus1.ito(Q_("m"))                                   # Centre fuselage y-coord
y_fus2 = ([0]*16)*Q_("m")
z_fus1 = [595, 714., 783, 812, 836, 843, 826, 801, 775, 749,
          724, 698, 672, 647, 621, 600] * Q_("mm")
z_fus1.ito(Q_("m"))
z_fus3 = [395, 115, 6.0, 9, 25, 41, 57, 73, 90, 106,
          122, 138, 154, 171, 187, 200] * Q_("mm")
z_fus3.ito(Q_("m"))
z_fus2 = (z_fus1 + z_fus3) / 2
Afusi = abs(y_fus1) * (z_fus1-z_fus2) * np.pi
SumAfus = sum(Afusi)
W_fus = Geometry.Masses.W_fus
W_ifus = W_fus * Afusi/SumAfus

rho_12 = 1
rho_23 = 1
r2 = rho_12 * np.sqrt((z_fus1-z_fus2)**2+(y_fus2-y_fus1)**2)
r4 = rho_23 * np.sqrt((z_fus1 - z_fus2)**2 + (y_fus2 - y_fus1)**2)
r = np.array([[z_fus1.magnitude - z_fus2.magnitude],
              [r2.magnitude],
              [y_fus2.magnitude - y_fus1.magnitude],
              [r4.magnitude],
              [z_fus2.magnitude-z_fus3.magnitude],
              [r4.magnitude],
              [y_fus2.magnitude - y_fus1.magnitude],
              [r2.magnitude]])
r *= Q_("m")
r = r[:,0]
theta = np.array([[np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [m.radians(90) - np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [m.radians(90) - np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [m.radians(90) - np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [m.radians(90) - np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))]])


alpha = np.array([[theta[0,0]],
                 [(theta[0,0]+theta[1,0])/2],
                 [(theta[0,0]+theta[1,0])/2],
                 [(theta[2,0] + theta[3,0])/2],
                 [theta[3,0]],
                 [(theta[2,0] + theta[3,0])/2],
                 [(theta[0,0]+theta[1,0])/2],
                 [(theta[0,0]+theta[1,0])/2]])
alpha = alpha[:,0]
ycg = np.array([[y_fus1],
                [r[1]*np.cos(theta[1])],
                [y_fus2 - y_fus1],
                [r[3]*np.sin(theta[3])],
                [y_fus1],
                [-(r[3]*np.sin(theta[3]))],
                [-(y_fus2 - y_fus1)],
                [-(r[1]*np.cos(theta[1]))]])
ycg = ycg[:,0]
zcg = np.array([[z_fus1],
                [z_fus2 + r[1]*np.sin(theta[1])],
                [z_fus2],
                [z_fus2 - r[3]*np.cos(theta[3])],
                [z_fus2],
                [z_fus2 - r[3]*np.cos(theta[3])],
                [z_fus2 + r[1]*np.sin(theta[1])]])
zcg = zcg[:,0]

s_pm = r * alpha
mpm = W_ifus * s_pm/(sum(s_pm))
