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
alpha_0 = 0                 # update
Z_cg = 0                    # Update
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
X_w = Geometry.CG.XLEMAC + 0.25 * Geometry.Wing.MAC      # CoP of main wing
X_cg = Geometry.CG.CG_mtow
mu_c = MTOW/(rho_a*Cbar*S_wing)
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
# specific parameters (only check if either Lambda or A changes)


# Iteration values
S_h = Geometry.H_tail.S
l_h = Geometry.H_tail.X_h - X_cg

# Values from graphs update with changing A and taper!!!
C_l_beta_CL_A = -0.02
dCnp_cl = -0.05
dCnp_CD0 = 11
CYr_cl2 = 0
ybar_r_b = 0.5
ybar_n_b = 0.55

# Calculated value
C_L = (MTOW*g0)/(0.5*rho_a*V_a**2*S_wing)


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
Re_sp = -B_sp/(2*A_sp)
Img_sp = m.sqrt(abs(Det_sp))/(2*A_sp)
Eigen_abs_sp = m.sqrt(Re_sp**2+Img_sp**2)
omega_0 = Eigen_abs_sp * (V_a/Cbar)
n_alpha = (CL_alpha * 0.5 * rho_a * V_a**2 * S_wing)/(MTOW * g0)
CAP = (omega_0**2)/n_alpha


# Phugoid

A_ph = 2 * mu_c * (CZ_alpha*CZq-2 * mu_c * Cm_alpha)
B_ph = 2 * mu_c * (CXu * Cm_alpha - Cmu * CX_alpha) + Cmq * \
       (CZu * CX_alpha - CXu * CZ_alpha)
C_ph = - CZ0 * (Cmu * CZ_alpha - CZu * Cm_alpha)
Det_ph = B_ph**2 - 4 * A_ph * B_ph
Re_ph = -B_ph/(2 * A_ph)
Img_ph = m.sqrt(abs(Det_ph))/(2*A_ph)
Damping_ph = -Re_ph/np.sqrt(Re_ph**2+Img_ph**2)
