# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

Stability derivatives range for level 1 flying qualities
"""

import sys
import math as m
import numpy as np
sys.path.append('../')

from Misc import Q_, ureg
from Geometry import Geometry
from Aerodynamics import Wing as Aero_wing
from Aerodynamics import HT as Aero_HT
from Aerodynamics import VT as Aero_VT
from Aerodynamics import General as Aero_gen
from Inertia import Inertia
from Performance import Performance

# This file will calculate the ranges for l_h,
# Xlemac, S_h and S_v in which StefX will havel level 1 flying qualities

# Input parameters
gamma_0 = 0                 # Assuming level flight
A = Geometry.Wing.A
             
Z_cg = Geometry.CG.Z_cg
CL_alpha = Aero_wing.CL_alpha
MTOW = Geometry.Masses.W_MTOW
S_wing = Geometry.Wing.S
taper = Geometry.Wing.taper
b = Geometry.Wing.b
V_a = (Performance.V_a_clean).magnitude
V_a *= Q_("1 m/s")
rho_a = Performance.rho_a.magnitude
rho_a *= Q_("1 kg/m**3")
g0 = Performance.g0.magnitude
g0 *= Q_("1 m/s**2 ")
Oswald_e = Aero_wing.Oswald_e
CNH_alpha = Aero_HT.C_Nh_alpha
dE_dalpha = Aero_wing.de_da
Vh_V = Aero_HT.Vh_v
Cbar = Geometry.Wing.MAC
I_yy = Inertia.I_yy
K_yy = np.sqrt(I_yy/MTOW)/(b)
CNW_alpha = Aero_wing.C_Nw_alpha
mu_c = MTOW/(rho_a*Cbar*S_wing)
mu_b = MTOW/(rho_a*b*S_wing)
CY_v_alpha = Aero_VT.C_Yv_alpha
dSig_dalpha = Aero_VT.q_dsigma_dbeta
Vv_V = Aero_VT.Vv_V
S_v = Geometry.V_tail.S
Z_v = Geometry.V_tail.Z_v
X_v = Geometry.V_tail.X_v
CD0 = Aero_gen.CD_0
I_xx = Inertia.I_xx
I_zz = Inertia.I_zz
K_xx = np.sqrt(I_xx/MTOW)/(b)
K_zz = np.sqrt(I_zz/MTOW)/(b)
CD0_alpha = 0                   # Assumed

# Values from graphs update with changing A and taper!!!
C_l_beta_CL_A = -0.02
dCnp_cl = -0.05
dCnp_CD0 = 11
CYr_cl2 = 0
ybar_r_b = 0.5
ybar_n_b = 0.55

# Calculated value
C_L = (MTOW*g0)/(0.5*rho_a*V_a**2*S_wing)
alpha_0 = C_L / CL_alpha


# CG value import

W_wing = Geometry.Masses.W_wing
W_htail = Geometry.Masses.W_htail
W_vtail = Geometry.Masses.W_vtail
W_fus = Geometry.Masses.W_fus
W_gear = Geometry.Masses.W_gear
W_engine = Geometry.Masses.W_engine
W_prop = Geometry.Masses.W_prop
W_fuelsys = Geometry.Masses.W_fuelsys
W_hydraulic = Geometry.Masses.W_hydraulic
W_flightcon = Geometry.Masses.W_flightcontrol
W_avionics = Geometry.Masses.W_avionics
W_elecsys = Geometry.Masses.W_elecsys
W_lehld = Geometry.Masses.W_lehld
W_flaperons = Geometry.Masses.W_flaperons
W_pilot = Geometry.Masses.W_pilot
W_fuel = Geometry.Masses.W_fuel
CG_wing_mac = Geometry.CG.CG_wing_mac
CG_fus = Geometry.CG.CG_fus
CG_lgear = Geometry.CG.CG_lgear
CG_engine = Geometry.CG.CG_engine
CG_prop = Geometry.CG.CG_prop
CG_fuelsys = Geometry.CG.CG_fuelsys
CG_hydraulics = Geometry.CG.CG_hydraulics
CG_flightcon = Geometry.CG.CG_flightcon
CG_avionics = Geometry.CG.CG_avionics
CG_elecsys = Geometry.CG.CG_elecsys
CG_lehld = Geometry.CG.CG_lehld
CG_flaperons = Geometry.CG.CG_flaperons
CG_pilot = Geometry.CG.CG_pilot
CG_fuel = Geometry.CG.CG_fuel
X_cg = Q_("1 m")                    # placeholders
X_w = Q_("1 m")
S_h = Q_("1 m**2")
l_h = Q_("5 m")
S_v = Q_("1 m**2")
# Local CG calculation for iterations

# Stability Derivatives (Longitudinal)

CX0 = (MTOW * g0)/(0.5 * rho_a * V_a**2 * S_wing) * m.sin(gamma_0)
CXu = -2 * C_L * m.tan(gamma_0)
CX_alpha = C_L * (1-(2*CL_alpha))/(Oswald_e*A*m.pi)
CZ0 = -(MTOW * g0)/(0.5*rho_a*V_a**2*S_wing)*m.cos(gamma_0)
CZu = -2*C_L
CZ_alpha = -CL_alpha
CZ_alphadot = -CNH_alpha * Vh_V**2 * dE_dalpha * S_h * l_h / (S_wing * Cbar)
CZq = -2 * Vh_V**2 * S_h * l_h / (S_wing * Cbar)
Cmu = 0
Cm_alpha = CNW_alpha * (X_cg - X_w) / Cbar - CNH_alpha * (1-dE_dalpha)*Vh_V**2\
            * S_h * l_h / (S_wing * Cbar)
Cm_alphadot = - CNH_alpha*(Vh_V)**2*S_h*l_h**2/(S_wing*Cbar**2)*dE_dalpha
Cmq = -1.1 * CNH_alpha * Vh_V**2 * (S_h*l_h**2)/(S_wing * Cbar**2)

# Stability Derivatives (Lateral)
z_arm = (((Z_v - Z_cg)/b)*np.cos(alpha_0)-((X_v-X_cg)/b)*np.sin(alpha_0))
x_arm = (((Z_v - Z_cg)/b)*np.sin(alpha_0)-((X_v-X_cg)/b)*np.cos(alpha_0)) #Check this
CYbeta_v = - CY_v_alpha * dSig_dalpha * Vv_V**2 * (S_v/S_wing)
Clbeta_w = C_L * C_l_beta_CL_A
Clbeta_v = CYbeta_v * z_arm
Cnbeta_w = C_L**2/(4*m.pi*A)
Cnbeta_v = CYbeta_v * x_arm
Clp_w = -((CL_alpha + CD0)*Cbar*b)/(24*S_wing)*(1+3*taper)
Cnp_w = dCnp_cl * C_L + dCnp_CD0 * CD0_alpha
Cnp_v = -2 * CYbeta_v * x_arm
CYr_v = 2 * CY_v_alpha * Vv_V**2 * (S_v*(X_v-X_cg))/(S_wing * b)
Clr_w = ybar_r_b**2 * C_L
Clr_v = CYr_v * z_arm
Cnr_w = ybar_n_b * C_L**2/(A*m.pi)
Cnr_v = CYr_v * x_arm

CYbeta = CYbeta_v
Clbeta = Clbeta_w + Clbeta_v
Cnbeta = Cnbeta_w + Cnbeta_v
Clp = Clp_w
Cnp = Cnp_w + Cnp_v
CYr = CYr_v
Clr = Clr_w + Clr_v
Cnr = Cnr_w + Cnr_v


# CAP

A_sp = -2 * mu_c * K_yy * (CZ_alphadot - 2 * mu_c)          # A for the Short period
B_sp = -CX_alpha * 2 * mu_c * K_yy + Cmq * (CZ_alphadot -2 * mu_c)\
        - CZq * Cm_alphadot - 2 * mu_c * Cm_alpha
C_sp = CZ_alpha * Cmq - CZq * Cm_alpha - 2 * mu_c * Cm_alpha
Det_sp = B_sp**2 - 4 * A_sp * B_sp
Re_sp = -B_sp/(2*A_sp)                          # Real part of eigenvalue SP
Img_sp = m.sqrt(abs(Det_sp))/(2*A_sp)           # Img part of eigenvalue SP
Eigen_abs_sp = m.sqrt(Re_sp**2+Img_sp**2)       # Abs value of eigenvalue SP
omega_0 = Eigen_abs_sp * (V_a/Cbar)             # Undamped nat. freq. SP
n_alpha = (CL_alpha * 0.5 * rho_a * V_a**2 * S_wing)/(MTOW * g0)    # loadfactor per AoA
CAP = (omega_0**2)/n_alpha


# Phugoid

A_ph = 2 * mu_c * (CZ_alpha*CZq-2 * mu_c * Cm_alpha)
B_ph = 2 * mu_c * (CXu * Cm_alpha - Cmu * CX_alpha) + Cmq * \
       (CZu * CX_alpha - CXu * CZ_alpha)
C_ph = - CZ0 * (Cmu * CZ_alpha - CZu * Cm_alpha)
Det_ph = B_ph**2 - 4 * A_ph * B_ph
Re_ph = -B_ph/(2 * A_ph)                        # Real part of eigenvalue Phug
Img_ph = m.sqrt(abs(Det_ph))/(2*A_ph)           # Img part of eigenvalue Phug
Damping_ph = -Re_ph/np.sqrt(Re_ph**2+Img_ph**2)  # Damping ratio Phug

# Aperiodic roll

A_ar = -4 * mu_b * K_xx
B_ar = Clp
Re_ar = Clp/(4 * mu_b * K_xx)
T_car = -1 / Re_ar * (b / V_a)

# Dutch roll (only stable, not level 1)

A_dr = 8 * mu_b**2 * K_zz
B_dr = - 2 * mu_b * (Cnr + 2 * K_zz * CYbeta)
RE_dr = -B_dr/(2 * A_dr)


