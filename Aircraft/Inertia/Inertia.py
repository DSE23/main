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
from Propulsion_and_systems import Engine

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
Y_mac = (b/6)*((1+2*taper_w)/(1+taper_w))           # Y location MAC
W_Htail = Geometry.Masses.W_htail
W_Vtail = Geometry.Masses.W_vtail
W_wing = Geometry.Masses.W_wing
W_wing = W_wing + Geometry.Masses.W_flaperons + Geometry.Masses.W_lehld
W_spar1 = StrucVal.Weightspar1
W_spar2 = StrucVal.Weightspar2
rho_rib = StrucVal.Density                          # Density rib material
k_rib = 0.5 * 10**-3                                # Value from Inertia report
t_ref = Q_("1 m")                                   # Reference thickness
t_c = Geometry.Wing.T_Cmax
t_rootwrib = c_rw * t_c                           # Thickness of airfoil at root
t_tipwrib = c_tw * t_c                            # Thickness of airfoil at tip
W_rib = k_rib * rho_rib * S_wing * (t_ref + (t_rootwrib + t_tipwrib)/2)
W_skin = W_wing - (W_spar1 + W_spar2 + W_rib)     # Skin weight
F_fs = W_spar1/W_wing                             # Front spar W fraction
F_rs = W_spar2/W_wing                             # Rear spar W fraction
F_skin = W_skin/W_wing                            # Skin W fraction
F_ribs = W_rib/W_wing                             # Ribs W fraction

# Fuselage
Z_cg = Geometry.CG.ZCG_mtow
Y_cg = 0                # Symmetry :P
X_cg = Geometry.CG.CG_mtow
x_fus = np.linspace(0, 6.2, 16) * Q_("m")
y_fus1 = [-400, -502, -520, -520, -509, -482, -439, -390, -342, -292,
          -243, -194, -145, -96, -47, -7.3] * Q_("-1 mm")  # left fuselage y-coord
y_fus1.ito(Q_("m"))                                   # Centre fuselage y-coord
y_fus2 = ([0]*16)*Q_("m")
Z_fus_origin = Geometry.CG.Z_fusorig
z_fus1 = [595, 714., 783, 812, 836, 843, 826, 801, 775, 749,
          724, 698, 672, 647, 621, 600] * Q_("mm")  # Fuselage Z-coord upper
z_fus1.ito(Q_("m"))
z_fus3 = [395, 115, 6.0, 9, 25, 41, 57, 73, 90, 106,
          122, 138, 154, 171, 187, 200] * Q_("mm")  # Fuselage Z-coord lower
z_fus3.ito(Q_("m"))
z_fus1 = Z_fus_origin - z_fus1
z_fus3 = Z_fus_origin - z_fus3
z_fus2 = (z_fus1 + z_fus3) / 2
#z_fus2 = ([0]*16)*Q_("m")
Afusi = abs(y_fus1) * (z_fus1-z_fus2) * np.pi      # Areas of fuselage sections
SumAfus = sum(Afusi)                                
W_fus = Geometry.Masses.W_fus
W_ifus = W_fus * Afusi/SumAfus


theta = np.array([[np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [Q_("90 deg") - np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [Q_("90 deg") - np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [Q_("90 deg") - np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [np.arctan((z_fus2-z_fus3)/(y_fus2-y_fus1))],
                 [Q_("90 deg") - np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))],
                 [np.arctan((y_fus2-y_fus1)/(z_fus1-z_fus2))]])
theta = theta[:, 0, :]


r1 = z_fus1 - z_fus2
r3 = y_fus1 - y_fus2
r5= z_fus2 - z_fus3
r2 = (abs(r1)*r3)/(np.sqrt(r1**2*(np.sin(theta[1]))**2+r3**2*(np.cos(theta[1]))**2))
r4 = (r5*r3)/(np.sqrt(r3**2*(np.sin(theta[3]))**2+r5**2*(np.cos(theta[3]))**2))

r = np.array([[r1.magnitude],
              [r2.magnitude],
              [r3.magnitude],
              [r4.magnitude],
              [r5.magnitude],
              [r4.magnitude],
              [r3.magnitude],
              [r2.magnitude]])
r *= Q_("m")
r = np.abs(r[:, 0])
alpha = np.array([[theta[0]],
                 [(theta[0] + theta[1])/2],
                 [(theta[0] + theta[1])/2],
                 [(theta[2] + theta[3])/2],
                 [theta[3]],
                 [(theta[2] + theta[3])/2],
                 [(theta[0] + theta[1])/2],
                 [(theta[0] + theta[1])/2]])
alpha = alpha[:, 0]
#ycgf = np.array([[y_fus1.magnitude],
#                [r[1].magnitude*np.cos(theta[1])],
#                [y_fus2.magnitude - y_fus1.magnitude],
#                [r[3].magnitude*np.sin(theta[3])],
#                [y_fus1.magnitude],
#                [-(r[3].magnitude*np.sin(theta[3]))],
#                [-(y_fus2.magnitude - y_fus1.magnitude)],
#                [-(r[1].magnitude*np.cos(theta[1]))]])
ycgf = np.array([[y_fus2.magnitude],
                [r[1].magnitude*np.cos(theta[1])],
                [y_fus2.magnitude - y_fus1.magnitude],
                [r[3].magnitude*np.sin(theta[3])],
                [y_fus2.magnitude],
                [-(r[3].magnitude*np.sin(theta[3]))],
                [-(y_fus2.magnitude - y_fus1.magnitude)],
                [-(r[1].magnitude*np.cos(theta[1]))]])

ycgf = ycgf[:, 0, :] * Q_("m")
zcgf = np.array([[z_fus1.magnitude],
                [z_fus2.magnitude - r[1].magnitude*np.sin(theta[1])],
                [z_fus2.magnitude],
                [z_fus2.magnitude + r[3].magnitude*np.cos(theta[3])],
                [z_fus3.magnitude],
                [z_fus2.magnitude + r[3].magnitude*np.cos(theta[3])],
                [z_fus2.magnitude],
                [z_fus2.magnitude - r[1].magnitude*np.sin(theta[1])]])

zcgf = zcgf[:, 0, :] * Q_("m")

xcgf = np.tile(x_fus, (8, 1)) * Q_("m")
s_pm = r * alpha
mpm = W_ifus * s_pm/(sum(s_pm))

zcgf_complete = (np.sum(np.sum(zcgf*mpm)))/(np.sum(np.sum(mpm)))
Z_CGF = Geometry.CG.ZCG_fus
if not abs(0.9*zcgf_complete) < abs(Z_CGF) < abs(1.1*zcgf_complete):
    print("Update ZCGf in Geometry to ", zcgf_complete)

I_xxpmf = mpm * ((ycgf - Y_cg)**2 + (zcgf - Z_cg)**2)     
I_yypmf = mpm * ((zcgf - Z_cg)**2 + (xcgf - X_cg)**2)
I_zzpmf = mpm * ((xcgf - X_cg)**2 + (ycgf - Y_cg)**2)
I_xzpmf = mpm * ((xcgf - X_cg) * (zcgf - Z_cg))
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
ycgw = np.tile(ycgw, (5, 1))
ycgw = ycgw * Q_("m")
ZCG_W = Geometry.CG.ZCG_wing
zcgw = np.ones(N_stw) * ZCG_W
a = -W_wing*(((c_rw*(1-taper_w))/sum(chordw))/N_stw)**2
C1 = 2/b.magnitude * (W_wing.magnitude/2 - b.magnitude**2/8 * a.magnitude)
C1 = C1 * Q_("m * kg")
A1 = b**2/(4*N_stw**2)*a/2 + b/(2*N_stw) * C1
B1 = (3 * b**2) / (2 * N_stw**2) * a/4 + b / (2 * N_stw) * C1
par = B1 - A1
A_airfoili = np.array([[(AreaAfoil(x_c0, x_c1, chordw)).magnitude],
                       [(AreaAfoil(x_c1, x_c2, chordw)).magnitude],
                       [(AreaAfoil(x_c2, x_c3, chordw)).magnitude],
                       [(AreaAfoil(x_c3, x_c4, chordw)).magnitude],
                       [(AreaAfoil(x_c4, x_c5-1*10**-10, chordw)).magnitude]])
A_airfoilfrac = A_airfoili/sum(A_airfoili)
N_windex = np.linspace(1, 40, 40)
mpmw = np.array([[((L1 * F_skin + A_airfoilfrac[0]*F_ribs)*(A1+par*(N_windex-1))).magnitude],
                 [(((L2 - L1) * F_skin + A_airfoilfrac[1] * F_ribs + F_fs)*(A1 + par*(N_windex-1))).magnitude],
                 [(((L3 - L2) * F_skin + A_airfoilfrac[2] * F_ribs)*(A1 + par*(N_windex-1))).magnitude],
                 [(((L4 - L3) * F_skin + A_airfoilfrac[3] * F_ribs + F_rs)*(A1 + par*(N_windex-1))).magnitude],
                 [(((1 - L4) * F_skin + A_airfoilfrac[4] * F_ribs)*(A1 + par* (N_windex-1))).magnitude]])
mpmw = mpmw[:, 0, 0, :] * Q_("kg")

I_xxpmw = mpmw * ((ycgw - Y_cg)**2 + (zcgw - Z_cg)**2)
I_yypmw = mpmw * ((zcgw - Z_cg)**2 + (xcgw - X_cg)**2)
I_zzpmw = mpmw * ((xcgw - X_cg)**2 + (ycgw - Y_cg)**2)
I_xzpmw = mpmw * ((xcgw - X_cg) * (zcgw - Z_cg))
I_xxw = 2 * sum(sum(I_xxpmw))
I_yyw = 2 * sum(sum(I_yypmw))
I_zzw = 2 * sum(sum(I_zzpmw))
I_xzw = 2 * sum(sum(I_xzpmw))


# Horizontal Tail
c_rh = Geometry.H_tail.c_r
b_h = Geometry.H_tail.b
c_th = Geometry.H_tail.c_t
Sweep_LE = Geometry.H_tail.Sweep_LE
C_ah = min(c_rh, ((b_h*np.tan(Sweep_LE))/2), (c_th + (b_h*np.tan(Sweep_LE))/2))
C_bh = b_h*np.tan(Sweep_LE)/2 + c_th
C_ch = max(c_rh, ((b_h*np.tan(Sweep_LE))/2), (c_th + (b_h*np.tan(Sweep_LE))/2))
rho_h = W_Htail / (0.5 * (-C_ah + C_bh + C_ch))
n_disctail = 25
x_h1 = np.linspace(0, C_ah, n_disctail)
x_h2 = np.linspace(C_ah, C_bh, n_disctail)
x_h3 = np.linspace(C_bh, C_ch, n_disctail)
x_h1 *= Q_("m")
hdx1 = x_h1[1]-x_h1[0]
x_h2 *= Q_("m")
hdx2 = x_h2[1]-x_h2[0]
x_h3 *= Q_("m")
hdx3 = x_h3[1]-x_h3[0]
y1 = rho_h/C_ah * x_h1
y2 = rho_h * (n_disctail*[1])
y3 = -((rho_h * x_h3)/(C_ch - C_bh)) + (rho_h * C_ch)/(C_ch - C_bh)
I_h = rho_h/12 * (-C_ah**3 + C_bh**3 + C_ch**2 * C_bh + C_ch*C_bh**2 + C_ch**3)
K_0 = 0.771                             # From source
sigmx_h = sum(y1 * hdx1 * x_h1) + sum(y2 * hdx2 * x_h2) + sum(y3 * hdx3 * x_h3)
sigm_h = sum(y1 * hdx1) + sum(y2 * hdx2) + sum(y3 * hdx3)
I_0yh = K_0 * (I_h - sigmx_h**2/sigm_h)
dxLE_h = 0.25 * c_rh - 0.25 * c_th
y_h = (2 * c_th * dxLE_h + c_th**2 + dxLE_h * c_rh +\
       c_th * c_rh + c_rh**2)/(3 * (c_rh + c_th))
H_rollcoeff = y_h/(b_h/6*((c_rh + 2 * c_th)/(c_rh + c_th)))
if not 0.95 < H_rollcoeff < 0.97:
    print(" !!!!Change k_4!!!!!")
k_4 = 0.88                              # From graphs, dependent on H_rollcoeff
I_0xh = (W_Htail * b_h**2 * k_4)/24 * ((c_rh+3 * c_th)/(c_rh + c_th))
I_0zh = I_0yh + I_0xh
ycgh = 0
zcgh = Geometry.CG.ZCG_htail             # Symmetric wing, so cg in middle of height
xcgh = Geometry.CG.CG_htail             # CG location in x-axis of the H-tail
I_xxh = I_0xh + W_Htail * ((ycgh - Y_cg)**2 + (zcgh - Z_cg)**2)
I_yyh = I_0yh + W_Htail * ((xcgh - X_cg)**2 + (zcgh - Z_cg)**2)
I_zzh = I_0zh + W_Htail * ((xcgh - X_cg)**2 + (ycgh - Y_cg)**2)
I_xzh = W_Htail * (xcgh - X_cg)*(zcgh - Z_cg)

# Vertical Tail
c_rv = Geometry.V_tail.c_r
b_v = Geometry.V_tail.b
c_tv = Geometry.V_tail.c_t
Sweep_LEv = Geometry.V_tail.Sweep_LE
C_av = min(c_rv, ((b_v*np.tan(Sweep_LEv))/2), (c_tv + (b_v*np.tan(Sweep_LEv))/2))
C_bv = c_tv + (b_v*np.tan(Sweep_LEv))/2
C_cv = max(c_rv, ((b_v*np.tan(Sweep_LEv))/2), (c_tv + (b_v*np.tan(Sweep_LEv))/2))
rho_v = W_Vtail/(0.5 * (-C_av + C_bv + C_cv))
x_v1 = np.linspace(0, C_av, n_disctail)
x_v2 = np.linspace(C_av, C_bv, n_disctail)
x_v3 = np.linspace(C_bv, C_cv, n_disctail)
x_v1 *= Q_("m")
vdx1 = x_v1[1] - x_v1[0]
x_v2 *= Q_("m")
vdx2 = x_v2[1] - x_v2[0]
x_v3 *= Q_("m")
vdx3 = x_v3[1] - x_v3[0]
yv1 = rho_v/C_av * x_v1
yv2 = rho_v * (n_disctail*[1])
yv3 = -((rho_v * x_v3)/(C_cv - C_bv)) + (rho_v * C_cv)/(C_cv - C_bv)
I_v = rho_v/12 * (-C_av**3 + C_bv**3 + C_cv**2 * C_bv + C_cv * C_bv**2 + C_cv**3)
dxLE_v = 0.25 * c_rv - 0.25*c_th
z_vbar = (2*c_tv* dxLE_v + c_tv**2 + dxLE_v * c_rv +\
       c_tv * c_rv + c_rv**2)/(3*(c_rv + c_tv)) 
V_rollcoef = z_vbar/((b_v/3)*(c_rv + 2* c_tv)/(c_rv + c_tv))
k_5 = 1.4
if not 1.6 < V_rollcoef < 1.65:
    print(" !!!! Change K_5!!!!!")
I_0xv = (W_Vtail * b_v**2 * k_5)/18 * (1+ (2 * c_rv * c_tv)/(c_rv + c_tv)**2)
sigmx_v = sum(yv1 * vdx1 * x_v1) + sum(yv2 * vdx2 * x_v2) + sum(yv3 * vdx3 * x_v3)
sigm_v = sum(yv1 * vdx1) + sum(yv2 * vdx2) + sum(yv3 * vdx3)
I_0zv = K_0 * (I_v - sigmx_v**2/sigm_v)
I_0yv = I_0zv + I_0xv
ycgv = 0
zcgv = Geometry.CG.ZCG_vtail
xcgv = Geometry.CG.CG_vtail
I_xxv = I_0xv + W_Vtail * ((ycgv - Y_cg)**2 + (zcgv - Z_cg)**2)
I_yyv = I_0yv + W_Vtail * ((xcgv - X_cg)**2 + (zcgv - Z_cg)**2)
I_zzv = I_0zv + W_Vtail * ((xcgv - X_cg)**2 + (ycgv - Y_cg)**2)
I_xzv = W_Vtail * (xcgv - X_cg)*(zcgv - Z_cg)

# Engine

W_engine = Engine.mass
I_0xe = Engine.ixg
I_0ye = Engine.iyg
I_0ze = Engine.izg
zcge = Geometry.CG.ZCG_engine
ycge = Q_("0 m")
xcge = Geometry.CG.CG_engine
I_xxe = I_0xe + W_engine * ((ycge - Y_cg)**2 + (zcge - Z_cg)**2)
I_yye = I_0ye + W_engine * ((xcge - X_cg)**2 + (zcge - Z_cg)**2)
I_zze = I_0ze + W_engine * ((xcge - X_cg)**2 + (ycge - Y_cg)**2)
I_xze = W_engine * (xcge - X_cg) * (zcge - Z_cg)


# Fuel

W_fuel = Geometry.Masses.W_fuel
xcgfuel = Geometry.CG.CG_fuel
zcgfuel = Geometry.CG.ZCG_fuel
ycgfuel = Q_("0 m")
I_xxfuel = W_fuel * ((ycgfuel - Y_cg)**2 + (zcgfuel - Z_cg)**2)
I_yyfuel = W_fuel * ((xcgfuel - X_cg)**2 + (zcgfuel - Z_cg)**2)
I_zzfuel = W_fuel * ((xcgfuel - X_cg)**2 + (ycgfuel - Y_cg)**2)
I_xzfuel = W_fuel * (xcgfuel - X_cg) * (zcgfuel - Z_cg)

# Pilot

W_pilot = Geometry.Masses.W_pilot
xcgp = Geometry.CG.CG_pilot
zcgp = Geometry.CG.ZCG_pilot
ycgp = Q_("0 m")
I_xxp = W_pilot * ((ycgp - Y_cg)**2 + (zcgp - Z_cg)**2)
I_yyp = W_pilot * ((xcgp - X_cg)**2 + (zcgp - Z_cg)**2)
I_zzp = W_pilot * ((xcgp - X_cg)**2 + (ycgp - Y_cg)**2)
I_xzp = W_pilot * (xcgp - X_cg) * (zcgp - Z_cg)

# Landing Gear

W_lg = Geometry.Masses.W_gear
xcglg = Geometry.CG.CG_lgear
ycglg = Q_("0 m")
zcglg = Geometry.CG.ZCG_lgear
I_xxlg = W_lg * ((ycglg - Y_cg)**2 + (zcglg - Z_cg)**2)
I_yylg = W_lg * ((xcglg - X_cg)**2 + (zcglg - Z_cg)**2)
I_zzlg = W_lg * ((xcglg - X_cg)**2 + (ycglg - Y_cg)**2)
I_xzlg = W_lg * (xcglg - X_cg) * (zcglg - Z_cg)

# Electronic Systems

W_elec = Geometry.Masses.W_elecsys
xcgel = Geometry.CG.CG_elecsys
zcgel = Geometry.CG.ZCG_elecsys
ycgel = Q_("0 m")
I_xxel = W_elec * ((ycgel - Y_cg)**2 + (zcgel - Z_cg)**2)
I_yyel = W_elec * ((xcgel - X_cg)**2 + (zcgel - Z_cg)**2)
I_zzel = W_elec * ((xcgel - X_cg)**2 + (ycgel - Y_cg)**2)
I_xzel = W_elec * (xcgel - X_cg) * (zcgel - Z_cg)

# Flight Controls

W_fcon = Geometry.Masses.W_flightcontrol
xcgc = Geometry.CG.CG_flightcon
zcgc = Geometry.CG.ZCG_flightcon
ycgc = Q_("0 m")
I_xxc = W_fcon * ((ycgc - Y_cg)**2 + (zcgc - Z_cg)**2)
I_yyc = W_fcon * ((xcgc - X_cg)**2 + (zcgc - Z_cg)**2)
I_zzc = W_fcon * ((xcgc - X_cg)**2 + (ycgc - Y_cg)**2)
I_xzc = W_fcon * (xcgc - X_cg) * (zcgc - Z_cg)


# Total new Inertia

I_xxnew = I_xxf + I_xxw + I_xxv + I_xxh + I_xxe + I_xxfuel + I_xxp + I_xxlg +\
          I_xxel + I_xxc
I_yynew = I_yyf + I_yyw + I_yyv + I_yyh + I_yye + I_yyfuel + I_yyp + I_yylg +\
          I_yyel + I_yyc
I_zznew = I_zzf + I_zzw + I_zzv + I_zzh + I_zze + I_zzfuel + I_zzp + I_zzlg +\
          I_zzel + I_zzc
I_xznew = I_xzf + I_xzw + I_xzv + I_xzh + I_xze + I_xzfuel + I_xzp + I_xzlg +\
          I_xzel + I_xzc
print(I_xxnew, "instead of", I_xx)
print(I_yynew, "instead of", I_yy)
print(I_zznew, "instead of", I_zz)
print(I_xznew, "instead of", 0)
I_xx = I_xxnew
I_yy = I_yynew
I_zz = I_zznew
I_xz = I_xznew
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(8):
    ax.scatter(xcgf[i], ycgf[i], zcgf[i])
for i in range(5):
    ax.scatter(xcgw[i], ycgw[i], zcgw)
    ax.scatter(xcgw[i], -ycgw[i], zcgw)
ax.set_xlim(0, 7)
ax.set_ylim(-2*(16/9),2*(16/9))
ax.set_zlim(-2,2)
#plt.scatter(x_fus, z_fus2)
#plt.scatter(x_fus, z_fus1)
#plt.scatter(x_fus, z_fus3)
plt.show()
