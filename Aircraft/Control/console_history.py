import sys, os
stdout_old = sys.stdout
sys.stdout = open(os.devnull, 'w')
sys.path.append('../')              # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
from Geometry import Geometry
from Aerodynamics import Wing as Awing
from Inertia import Inertia
from Performance import Performance
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import pandas as pd
import math as m
import time
sys.stdout = stdout_old
t0 = time.time()

np.seterr(all='raise')

# Variables
l_a = Q_("0.405 m")        # Set aileron length
l_e = Q_("0.2    m")        # Set elevator length
cr_c= Q_("0.2     ")
n_of_disc_w = 20            # number of parts wing is discretized
n_of_disc_h = 10            # number of parts HT is discretized
n_of_disc_v = 10            # number of parts VT is discretized
da = Q_("0 deg")            # aileron deflection
dr = Q_("0 deg")            # rudder deflection
de = Q_("0 deg")            # elevator deflection
alpha_nose = Q_("0.05 rad") # angle of attack of nose
beta_nose  = Q_("0. rad")   # angle of sideslip of nose
V_inf = Q_("118 m/s")     # V infinity
t_current = Q_("0.0 s")       # Start time of sim
dt = Q_("0.01 s")           # Time step of sim
t_end = Q_("30. s")         # End time of sim
l_h = Q_("3.6444 m")        # Tail arm ac-ac horizontal
l_v = Q_("3.7 m")           # Tail arm ac-ac vertical
p = Q_("0. 1/s")            # initial roll rate  [rad/s]
q = Q_("0. 1/s")            # initial pitch rate [rad/s]
r = Q_("0. 1/s")            # initial yaw rate   [rad/s]
Phi   = Q_("0. rad")        # Initial euler angle around x-axis
Psi   = Q_("0. rad")        # Initial euler angle around z-axis
Theta = Q_("0. rad")        # Initial euler angle around y-axis
lin_ran_alpha = 4          # Linear range of angle of attack and elevator defl.
w = Q_("0. m/s")
u = V_inf
v = Q_("0. m/s")
gamma = Q_("0. rad")
Xi = Q_("0 rad")
mtow = Geometry.Masses.W_MTOW
g0 = Performance.g0.magnitude * Q_("m/s**2")
W = mtow*g0

# Import Aircraft Geometry
I_yy = Inertia.I_yy
I_xx = Inertia.I_xx
I_zz = Inertia.I_zz
I_xz = Inertia.I_xz
I_star = I_xx*I_zz-I_xz**2

b_w = Geometry.Wing.b
S_w = Geometry.Wing.S
c_r_w = Geometry.Wing.c_r
c_t_w = Geometry.Wing.c_t
AR_w = Geometry.Wing.A
e_w = Awing.Oswald_e
e_h = e_w
e_v = e_w
b_h = Geometry.H_tail.b
S_h = Geometry.H_tail.S
c_r_h = Geometry.H_tail.c_r
c_t_h = Geometry.H_tail.c_t
AR_h = Geometry.H_tail.A
b_v = Geometry.V_tail.b
S_v = Geometry.V_tail.S
c_r_v = Geometry.V_tail.c_r
c_t_v = Geometry.V_tail.c_t
AR_v = Geometry.V_tail.A
t_c_v = Geometry.V_tail.t_c
Z_v = Geometry.V_tail.Z_v

xcg    = Geometry.CG.CG_mtow
xlemac = Geometry.CG.XLEMAC
MAC = Geometry.Wing.MAC
MAC_htail = Geometry.H_tail.MAC
X_h = Geometry.CG.X_htail
X_w = Geometry.CG.X_wing

sweep_LE = Geometry.Wing.Sweep_LE.magnitude*1

cabin_width = Geometry.Fuselage.cabin_w + 0.1 * ureg.m  # import cabin width
vt_width = t_c_v * c_r_v + 0.1 * ureg.m                 # import width of VT
rho = Performance.rho_0.magnitude * Q_("kg/m**3")                        # import density

bloc_w = (b_w-cabin_width) / n_of_disc_w  # span of each station wing
bloc_h = b_h / n_of_disc_h  # Span of each station HT
bloc_v = b_v / n_of_disc_v  # Span of each station VT
half_b_w = b_w / 2          # half span wing
half_b_h = b_h / 2          # half span HT
y_mac = half_b_w*(c_r_w-MAC)/(c_r_w-c_t_w)

V_s = Performance.V_stall_clean            # stall speed
V_a = Performance.V_a_clean                  # manoeuvring speed

# Definitions
def local_chord(z, c_r, c_t, half_b):
    # Calculates the chord at location z(distance from center)
    return c_r - (c_r - c_t) / half_b * z

data = pd.read_csv('aerodynamic_data_ms15.dat', ' ', header=None).values

def lookup_data(alpha, ca_c, da):
    # Looksup and interpolates Cl and Cd based on alpha, ca_c and da
    alpha = round(alpha,5)
    alpha = m.degrees(alpha)
    indexca_c = int(100*ca_c)-1
    if da == 0.:
        indexda = abs(da)//5
        localdata = data[int(indexda*50*51+indexca_c*51):int(indexda*50*51+(indexca_c+1)*51),:]
        non_zero_max = max(np.argwhere(localdata[:, 0]))[0]  # last non-zero row
        localdata = localdata[:non_zero_max+1,:]
    else:
        index1da = abs(da)//5
        index2da = abs(da)//5 + 1
        index3da = 0
        localdata1 = data[int(index1da*50*51+indexca_c*51):int(index1da*50*51+(indexca_c+1)*51),:]
        localdata2 = data[int(index2da*50*51+indexca_c*51):int(index2da*50*51+(indexca_c+1)*51),:]
        localdata3 = data[int(index3da*50*51+indexca_c*51):int(index3da*50*51+(indexca_c+1)*51),:]
        non_zero_max1 = max(np.argwhere(localdata1[:, 0]))[0]  # last non-zero row
        non_zero_max2 = max(np.argwhere(localdata2[:, 0]))[0]  # last non-zero row
        non_zero_max3 = max(np.argwhere(localdata3[:, 0]))[0]  # last non-zero row
        non_zero_max = min((non_zero_max1,non_zero_max2,non_zero_max3))
        
        localdata1 = localdata1[:non_zero_max+1,:]
        localdata2 = localdata2[:non_zero_max+1,:]    
        localdata3 = localdata3[:non_zero_max+1,:]   
        
        localdata_alpha = localdata2[:,0]
        localdata_cl    = (localdata2[:,1] - localdata1[:,1])/5 *abs(da) + localdata3[:,1]
        localdata_cd    = (localdata2[:,2] - localdata1[:,2])/5 *abs(da) + localdata3[:,2]
        localdata_cm    = (localdata2[:,4] - localdata1[:,4])/5 *abs(da) + localdata3[:,4]
    
        localdata = np.array([localdata_alpha,
                              localdata_cl,
                              localdata_cd,
                              localdata_cd,
                              localdata_cm])
        localdata = np.transpose(localdata)
    
    Cl_local = interpolate.interp1d(localdata[:,0], localdata[:,1], 'linear', fill_value='extrapolate')
    Cd_local = interpolate.interp1d(localdata[:,0], localdata[:,2], 'linear', fill_value='extrapolate')
    Cm_local = interpolate.interp1d(localdata[:,0], localdata[:,4], 'linear', fill_value='extrapolate')
    if da >= 0:
        Cl = Cl_local(alpha)
        Cd = Cd_local(alpha)
        Cm = Cm_local(alpha)
    else:
        Cl = -Cl_local(-alpha)
        Cd = Cd_local(-alpha)
        Cm = -Cm_local(-alpha)
    
    Cn_lookup = m.cos(m.radians(alpha))*Cl + m.sin(m.radians(alpha))*Cd
    if Cn_lookup !=0:
        xcp = 0.25-Cm/Cn_lookup
    else:
        xcp = 0.25
    return Cl, Cd, Cm, xcp

"============================================================================="

de = -5.9
da = 29.5

alpha_w_l = alpha_nose
alpha_w_r = alpha_nose
alpha_h = alpha_nose
V_local = V_inf * m.cos(beta_nose)

n_chords_w = int(n_of_disc_w / 2)
kwlst = n_chords_w * [0] + [-cabin_width / 2, 0 * ureg.m, cabin_width / 2] + n_chords_w * [0]
for j in range(n_chords_w):
    kwlst[j] = bloc_w * j - half_b_w
    kwlst[-1 - j] = half_b_w - bloc_w * j
disc_wing_w = np.zeros((len(kwlst)-1, 10))
disc_wing_w[(range(n_chords_w)), 0] = da  
disc_wing_w[(range(int(n_of_disc_w + 2 - n_chords_w), n_of_disc_w + 2)), 0] = -da  # set aileron for negative stations

plst = []
tlst = np.arange(0,2,dt.magnitude)

for t_current in tlst:
    #Wing
    for i in range(0, len(kwlst)-1):
        da_local = disc_wing_w[i,0]                        # Local aileron deflection of piece
        b1 = kwlst[i]                                       # Y boundary left of piece
        b2 = kwlst[i+1]                                     # Y boundary right of piece
        y_i = (b1 + b2) / 2                                 # Y centre of piece
        c1 = local_chord(abs(b1), c_r_w, c_t_w, half_b_w)   # Chord left of piece
        c2 = local_chord(abs(b2), c_r_w, c_t_w, half_b_w)   # Chord right of piece
        cc = local_chord(abs(y_i), c_r_w, c_t_w, half_b_w)  # Chord centre of piece
        ca_c = (l_a / cc).magnitude                         # percentage aileron chord over local aileron
        Sloc = (c1 + c2) / 2 * (b2 - b1)                    # Surface area of piece
        
        roll_induced_alpha = p * y_i / V_local
        alpha_w = alpha_nose + roll_induced_alpha 
        
        Cl, Cd, Cm, xcp = lookup_data(alpha_w, ca_c, da_local)
        Cn_w = -Cl*m.cos(alpha_w) - Cd*m.sin(alpha_w)
        Ct_w = Cl*m.sin(alpha_w) - Cd*m.cos(alpha_w) 
        
        Fn_w_d = 0.5 * rho * V_local ** 2 * Sloc * Cn_w
        Ft_w_d = 0.5 * rho * V_local ** 2 * Sloc * Ct_w
        
        disc_wing_w[i,1] = Fn_w_d.magnitude
        disc_wing_w[i,2] = Ft_w_d.magnitude
        disc_wing_w[i,3] = (Fn_w_d*y_i).magnitude
        

    # Wing simple
    
    dx_w = (y_mac+0.25*MAC)-xcg
    Cl, Cd, Cm, xcp = lookup_data(alpha_nose, 0.2,0.)
    Cn_w = -Cl*m.cos(alpha_nose) - Cd*m.sin(alpha_nose)
    Ct_w = Cl*m.sin(alpha_nose) - Cd*m.cos(alpha_nose) 
    
    Fn_w = 0.5 * rho * V_local ** 2 * S_w * Cn_w
    Ft_w = 0.5 * rho * V_local ** 2 * S_w * Ct_w

    #Horizontal Tail
    dx_h = l_h
    Cl, Cd, Cm, xcp = lookup_data(alpha_h, 0.5,de)
    Cn_h = -Cl*m.cos(alpha_h) - Cd*m.sin(alpha_h)
    Ct_h = Cl*m.sin(alpha_h) - Cd*m.cos(alpha_h) 
    
    Fn_h = 0.5 * rho * V_local ** 2 * S_h * Cn_h
    Ft_h = 0.5 * rho * V_local ** 2 * S_h * Ct_h
    
    # Moments and Forces
    sum_Fn_w = sum(disc_wing_w[:,1])*Q_("N")
    My = dx_w * sum_Fn_w + dx_h * Fn_h
    Mx = sum(disc_wing_w[:,3])*Q_("N*m")
    
    p_dot = Mx/I_xx
    q_dot = My/I_yy
    
    p += p_dot * dt
    q += q_dot * dt
    
    plst.append(m.degrees(p.magnitude))

plt.plot(tlst,plst)
print("Fz:",Fn_w+Fn_h+W)