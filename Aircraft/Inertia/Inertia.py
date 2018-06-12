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
from Structures import Wing as strucwing

I_yy = Q_("1594 kg m**2")
I_xx = Q_("1089 kg m**2")
I_zz = Q_("2629 kg m**2")


# New inertia Calculation
Z_w = 0                     # Check this!!!!!!!!!!
b = Geometry.Wing.b
c_rw = Geometry.Wing.c_r
taper_w = Geometry.Wing.taper
Sparloc1 = strucwing.ChSpar1
Sparloc2 = strucwing.ChSpar2
XLEMAC = Geometry.CG.XLEMAC
Sweep_LE = Geometry.Wing.Sweep_LE
Y_mac = (b/6)*((1+2*taper_w)/(1+taper))
# Fuselage
Z_cg = Geometry.CG.Z_cg
Y_cg = 0                #Symmetry :P
X_cg = Geometry.CG.CG_mtow
x_fus = np.linspace(0, 6.2, 16) * Q_("m")
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
r = r[:, 0]
theta = np.array([[np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [m.radians(90) - np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [m.radians(90) - np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [m.radians(90) - np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [m.radians(90) - np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))]])
theta = theta[:,0,:] * Q_("rad")

alpha = np.array([[theta[0]],
                 [(theta[0] + theta[1])/2],
                 [(theta[0] + theta[1])/2],
                 [(theta[2] + theta[3])/2],
                 [theta[3]],
                 [(theta[2] + theta[3])/2],
                 [(theta[0] + theta[1])/2],
                 [(theta[0] + theta[1])/2]])
alpha = alpha[:, 0]
ycgf = np.array([[y_fus1.magnitude],
                [r[1].magnitude*np.cos(theta[1])],
                [y_fus2.magnitude - y_fus1.magnitude],
                [r[3].magnitude*np.sin(theta[3])],
                [y_fus1.magnitude],
                [-(r[3].magnitude*np.sin(theta[3]))],
                [-(y_fus2.magnitude - y_fus1.magnitude)],
                [-(r[1].magnitude*np.cos(theta[1]))]])
ycgf = ycgf[:, 0, :] * Q_("m")
zcgf = np.array([[z_fus1.magnitude],
                [z_fus2.magnitude + r[1].magnitude*np.sin(theta[1, 0])],
                [z_fus2.magnitude],
                [z_fus2.magnitude - r[3].magnitude*np.cos(theta[3, 0])],
                [z_fus2.magnitude],
                [z_fus2.magnitude - r[3].magnitude*np.cos(theta[3, 0])],
                [z_fus2.magnitude],
                [z_fus2.magnitude + r[1].magnitude*np.sin(theta[1, 0])]])
zcgf = zcgf[:, 0, :] * Q_("m")
xcgf = np.tile(x_fus, (8,1)) * Q_("m")
s_pm = r * alpha
mpm = W_ifus * s_pm/(sum(s_pm))
I_xxpmf = mpm * ((ycgf - Y_cg)**2 + (zcgf - Z_cg)**2)
I_yypmf = mpm * ((zcgf - Z_cg)**2 + (xcgf - X_cg)**2)
I_zzpmf = mpm * ((xcgf - X_cg)**2 + (ycgf - Y_cg)**2)
I_xzpmf = mpm * ((xcgf - X_cg)+(zcgf - Z_cg))
I_xxf = np.sum(np.sum(I_xxpmf))
I_yyf = np.sum(np.sum(I_yypmf))
I_zzf = np.sum(np.sum(I_zzpmf))
I_xzf = np.sum(np.sum(I_xzpmf))

# Wing

def AreaAfoil(x1, x2, chord):
    n = 100
    dx = (x2 - x1)/n
    dxlength = dx * chord
    area_cell = 0
    x = x1
    for i in range(n):
        x = x + dx
        area_cell = dxlength* strucwing.airfoilordinate(x)
        area_cell = area_cell * 2
        print(area_cell)
    return area_cell

N_stw = 40
ycgw = []
chordw = []
xapexw = []
for i in range(N_stw):
    ycgw = np.append(ycgw, ((i+1)/N_stw * b/2 - 1/2 * (b/2)/N_stw))
    chordw = np.append(chordw, (c_rw * (1- (((i+1)-0.5) * (1-taper_w))/N_stw)))
    xapexw = np.append(xapexw, (XLEMAC + ((ycgw[i]*Q_("m")) - Y_mac) * np.tan(Sweep_LE)))
ycgw *= Q_("m")
chordw *= Q_("m")
xapexw *= Q_("m")
x_c0 = 0
x_c1 = 0.5 * Sparloc1
x_c2 = 0.5 * (Sparloc1 + Sparloc2)
x_c3 = 0.5 * (1 + Sparloc2)
x_c4 = 1.0
X_cg = np.array([[xapexw],
                [1],
                [1],
                [1],
                [1]])