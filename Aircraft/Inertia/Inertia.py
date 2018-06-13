# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

Inertia parameters

"""
import sys
sys.path.append('../')   # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import math as m
import numpy as np
from Geometry import Geometry
from Structures import Wing as strucwing
from Structures import StrucVal

I_yy = Q_("1594 kg m**2")
I_xx = Q_("1089 kg m**2")
I_zz = Q_("2629 kg m**2")


# New inertia Calculation
b = Geometry.Wing.b
S_wing = Geometry.Wing.S
c_rw = Geometry.Wing.c_r
taper_w = Geometry.Wing.taper
c_tw = c_rw * taper_w
Sparloc1 = strucwing.ChSpar1
Sparloc2 = strucwing.ChSpar2
XLEMAC = Geometry.CG.XLEMAC
Sweep_LE = Geometry.Wing.Sweep_LE
Y_mac = (b/6)*((1+2*taper_w)/(1+taper_w))
W_Htail = Geometry.Masses.W_htail
W_wing = Geometry.Masses.W_wing
W_spar1 = StrucVal.Weightspar1
W_spar2 = StrucVal.Weightspar2
rho_rib = StrucVal.Density                          # Density rib material
k_rib = 0.5 * 10**-3                                # Value from Inertia report
t_ref = Q_("1 m")                                   # Reference thickness
t_c = Geometry.Wing.T_Cmax
t_rootwrib = c_rw * t_c
t_tipwrib = c_tw * t_c
W_rib = k_rib * rho_rib * S_wing * (t_ref + (t_rootwrib + t_tipwrib)/2)
W_skin = W_wing - (W_spar1 + W_spar2 + W_rib)
F_fs = W_spar1/W_wing
F_rs = W_spar2/W_wing
F_skin = W_skin/W_wing
F_ribs = W_rib/W_wing

# Fuselage
Z_cg = Geometry.CG.Z_cg
Y_cg = 0                # Symmetry :P
X_cg = Geometry.CG.CG_mtow
x_fus = np.linspace(0, 6.2, 16) * Q_("m")
y_fus1 = [-400, -502, -520, -520, -509, -482, -439, -390, -342, -292,
          -243, -194, -145, -96, -47, -7.3] * Q_("mm")  # left fuselage y-coord
y_fus1.ito(Q_("m"))                                   # Centre fuselage y-coord
y_fus2 = ([0]*16)*Q_("m")
z_fus1 = [595, 714., 783, 812, 836, 843, 826, 801, 775, 749,
          724, 698, 672, 647, 621, 600] * Q_("mm")  # Fuselage Z-coord upper
z_fus1.ito(Q_("m"))
z_fus3 = [395, 115, 6.0, 9, 25, 41, 57, 73, 90, 106,
          122, 138, 154, 171, 187, 200] * Q_("mm")  # Fuselage Z-coord lower
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
xcgf = np.tile(x_fus, (8, 1)) * Q_("m")
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
        area_cell = dxlength * strucwing.airfoilordinate(x)
        area_cell = area_cell * 2
    return area_cell


N_stw = 40
ycgw = []
chordw = []
xapexw = []
for i in range(N_stw):
    ycgw = np.append(ycgw, ((i+1)/N_stw * b/2 - 1/2 * (b/2)/N_stw))
    chordw = np.append(chordw, (c_rw * (1 - (((i+1)-0.5) * (1-taper_w))/N_stw)))
    xapexw = np.append(xapexw, (XLEMAC + ((ycgw[i] * Q_("m")) - Y_mac) * np.tan(Sweep_LE)))
ycgw *= Q_("m")
chordw *= Q_("m")
xapexw *= Q_("m")
xms = (Sparloc1 + Sparloc2)/2
x_c0 = 0
x_c1 = 0.5 * Sparloc1
x_c2 = 0.5 * (Sparloc1 + xms)
x_c3 = 0.5 * (xms + Sparloc2)
x_c4 = 0.5 * (1 + Sparloc2)
x_c5 = 1.0
L1 = x_c1
L2 = x_c2 
L3 = x_c3 
L4 = x_c4 
L5 = 1
xcgw = np.array([[(xapexw + L1/2 * chordw).magnitude],
                [(xapexw + Sparloc1 * chordw).magnitude],
                [(xapexw + xms * chordw).magnitude],
                [(xapexw + Sparloc2 * chordw).magnitude],
                [(xapexw + (L4+1)/2 * chordw).magnitude]])
xcgw = xcgw[:, 0, :] * Q_("m")
ycgw = np.tile(ycgw, (5,1))
ycgw = ycgw * Q_("m")
zcgw = np.zeros(N_stw)
a = -W_wing*(((c_rw*(1-taper_w))/sum(chordw))/N_stw)**2
C1 = 2/b.magnitude *(W_wing.magnitude/2 - b.magnitude**2/8 * a.magnitude)
C1 = C1 * Q_("m * kg")
A1 = b**2/(4*N_stw**2)*a/2 + b/(2*N_stw) * C1
B1 = (3 * b**2)/(2*N_stw**2) * a/4 + b /(2 * N_stw) * C1
par = B1 - A1
A_airfoili = np.array([[(AreaAfoil(x_c0, x_c1, chordw)).magnitude],
                       [(AreaAfoil(x_c1, x_c2, chordw)).magnitude],
                       [(AreaAfoil(x_c2, x_c3, chordw)).magnitude],
                       [(AreaAfoil(x_c3, x_c4, chordw)).magnitude],
                       [(AreaAfoil(x_c4, x_c5-1*10**-10, chordw)).magnitude]])
# A_airfoili = A_airfoili[:, 0, :]*Q_("m**2")
A_airfoilfrac = A_airfoili/sum(A_airfoili)
# W_sec = ycgw/((b/2) * N_stw) * W_wing
# mpmw = W_sec * A_airfoilfrac
N_windex = np.linspace(1,40,40)
mpmw = np.array([[((L1 * F_skin + A_airfoilfrac[0]*F_ribs)*(A1+par*(N_windex-1))).magnitude],
                 [(((L2 - L1) * F_skin + A_airfoilfrac[1] * F_ribs + F_fs)*(A1 + par*(N_windex-1))).magnitude],
                 [(((L3 - L2) * F_skin + A_airfoilfrac[2] * F_ribs)*(A1 + par*(N_windex-1))).magnitude],
                 [(((L4 - L3) * F_skin + A_airfoilfrac[3] * F_ribs + F_rs)*(A1 + par*(N_windex-1))).magnitude],
                 [(((1 - L4) * F_skin + A_airfoilfrac[4] * F_ribs)*(A1 + par* (N_windex-1))).magnitude]])
mpmw = mpmw[:, 0, 0, :] * Q_("kg")

I_xxpmw = mpmw * ((ycgw - Y_cg)**2 + (zcgw - Z_cg)**2)
I_yypmw = mpmw * ((zcgw - Z_cg)**2 + (xcgw - X_cg)**2)
I_zzpmw = mpmw * ((xcgw - X_cg)**2 + (ycgw - Y_cg)**2)
I_xypmw = mpmw * ((xcgw - X_cg) + (zcgw - Z_cg))
I_xxw = 2 * sum(sum(I_xxpmw))
I_yyw = 2 * sum(sum(I_yypmw))
I_zzw = 2 * sum(sum(I_zzpmw))
I_xyw = 2 * sum(sum(I_xypmw))


# Horizontal Tail
c_rh = Geometry.H_tail.c_r
b_h = Geometry.H_tail.b
c_th = Geometry.H_tail.c_t
Sweep_LE = Geometry.H_tail.Sweep_LE
C_ah = min(c_rh, ((b_h*np.tan(Sweep_LE))/2), (c_th + (b_h*np.tan(Sweep_LE))/2) )
C_bh = b_h*np.tan(Sweep_LE)/2 + c_th
C_ch = max(c_rh, ((b_h*np.tan(Sweep_LE))/2), (c_th + (b_h*np.tan(Sweep_LE))/2) )
rho_h = W_Htail/ (0.5 * ( -C_ah + C_bh + C_ch))
x_h = np.linspace(0, c_rh, 50)
x_h1 = np.linspace(0, C_ah, 20)
x_h2 = np.linspace(C_ah, C_bh, 20)
x_h3 = np.linspace(C_bh, C_ch, 20)
x_h1 *= Q_("m")
hdx1 = x_h1[1]-x_h1[0]
x_h2 *= Q_("m")
hdx2 = x_h2[1]-x_h2[0]
x_h3 *= Q_("m")
hdx3 = x_h3[1]-x_h3[0]
y1 = rho_h/C_ah * x_h1
y2 = rho_h * (20*[1])
y3 = -((rho_h* x_h3)/(C_ch - C_bh)) + (rho_h * C_ch)/(C_ch - C_bh)
I_h = rho_h/12 * (-C_ah**3 + C_bh**3 + C_ch**2 * C_bh + C_ch*C_bh**2 + C_ch**3)
K_0 = 0.771                             # From source 
sigmx_h = sum(y1 * hdx1 * x_h1) + sum(y2 * hdx2 * x_h2) + sum(y3 * hdx3 * x_h3)
sigm_h = sum(y1 * hdx1) + sum(y2 * hdx2) + sum(y3 * hdx3)
I_0yh = K_0 * (I_h - sigmx_h**2/sigm_h)
y_h = Q_("0.6 m")
H_rollcoeff = y_h/(b_h/6*((c_rh + 2 * c_th)/(c_rh + c_th)))
k_4 = 0.88                              # From graphs, dependent on H_rollcoeff
I_0xh = (W_Htail * b_h**2 * k_4)/24 * ((c_rh+3*c_th)/(c_rh +c_th))
I_0zh = I_0yh + I_0xh
ycgh = 0
zcgh = Geometry.H_tail.Z_h              # Symmetric wing, so cg in middle of height
xcgh = Geometry.CG.CG_htail             # CG location in x-axis of the H-tail
I_xxh = I_0xh + W_Htail * ((ycgh - Y_cg)**2 + (zcgh - Z_cg)**2)
I_yyh = I_0yh + W_Htail * ((xcgh - X_cg)**2 + (zcgh - Z_cg)**2)
I_zzh = I_0zh + W_Htail * ((xcgh - X_cg)**2 + (ycgh - Y_cg)**2)
Ixz = W_Htail * (xcgh - X_cg)*(zcgh - Z_cg)

# Vertical Tail
c_rv = Geometry.V_tail.c_r
b_v = Geometry.V_tail.b
c_th = Geometry.V_tail.c_t
Sweep_LEv = Geometry.V_tail.Sweep_LE