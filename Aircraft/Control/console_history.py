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
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import pandas as pd
import math as m
import time
sys.stdout = stdout_old
t0 = time.time()

np.seterr(all='raise')

# Variables
l_a = Q_("0.405 m")        # Set aileron length
cr_c= Q_("0.45     ")
ce_c= Q_("0.5 ")
n_of_disc_w = 20            # number of parts wing is discretized
n_of_disc_h = 10            # number of parts HT is discretized
n_of_disc_v = 10            # number of parts VT is discretized
da = Q_("0 deg")            # aileron deflection
dr = Q_("0 deg")            # rudder deflection
de = Q_("0 deg")            # elevator deflection
alpha_nose = Q_("0.012 rad") # angle of attack of nose
beta_nose  = Q_("0. rad")   # angle of sideslip of nose
V_inf = Q_("50 m/s")     # V infinity
t_current = Q_("0.0 s")       # Start time of sim
dt = Q_("0.001 s")           # Time step of sim
t_end = Q_("30. s")         # End time of sim
p = Q_("0. 1/s")            # initial roll rate  [rad/s]
q = Q_("0. 1/s")            # initial pitch rate [rad/s]
r = Q_("0. 1/s")            # initial yaw rate   [rad/s]
Phi   = Q_("0. rad")        # Initial euler angle around x-axis
Psi   = Q_("0. rad")        # Initial euler angle around z-axis
Theta = Q_("0. rad")        # Initial euler angle around y-axis
lin_ran_alpha = Q_("10 deg")          # Linear range of angle of attack and elevator defl.
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

X_h = Geometry.H_tail.X_h
X_v = Geometry.V_tail.X_v
MAC_h = Geometry.H_tail.MAC
MAC_v = Geometry.V_tail.MAC
l_h = X_h + 0.25 * MAC_h
l_v = X_v + 0.25 * MAC_v

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
    alpha = alpha.to('degree')
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


V = Q_("99 m/s")
def trimlift(angles):
    alpha = angles[0]
    de = angles[1]
    Cl, Cd, Cm, xcp = lookup_data(alpha, 0.3,0.)
    Clh, Cdh, Cmh, xcph = lookup_data(alpha, 0.5, de)
    Cn_w = -Cl*np.cos(alpha_w) - Cd*np.sin(alpha_w)
    Cn_h = -Clh*np.cos(alpha_w) - Cd*np.sin(alpha_w)
    q_dp = 0.5*rho*V**2
    L_w = - Cn_w * q_dp * S_w
    L_h = - Cn_h * q_dp * S_h
    F_z = W  - L_w - L_h
    x_w = xlemac + 0.25 * MAC
    x_h = X_h + 0.25 * MAC_htail
    M_y = -L_w * (x_w - xcg) + Cm * q_dp * S_w * MAC + Cmh * q_dp * S_h * MAC_htail +\
          -L_h * (x_h - xcg)
    return M_y.magnitude**2 + F_z.magnitude**2      

def trim():
    
    initial_guess = [0.02, -2]
    result = optimize.minimize(trimlift, initial_guess, method='POWELL', options={'ftol': 0.00001})
    print (result.success)
    angles =  result.x
    alpha_t = angles[0] * Q_("rad")
    de_t = angles[1] * Q_("deg")
    Cltrimw, Cdtrimw, Cmtrimw, xcp = lookup_data(alpha_t, ca_c, 0)
    Cltrimh, Cdtrimh, Cmtrimh, xcp = lookup_data(alpha_t, ce_c, de_t)
    q_dp = 0.5 * rho * V**2
    Treq = (np.cos(alpha_t)*(Cdtrimw*S_w + Cdtrimh * S_h)+np.sin(alpha_t)*
            (Cltrimw * S_w + Cltrimh * S_h)) * q_dp
    
    return alpha_t, de_t, Treq   
    
#def trim(V):
#    lin_ran_de = Q_("10 deg")
#    alpha_rad = 0.087
#    Claw1, Cdaw0, Cmaw1, xcp = lookup_data(Q_("0 deg"), ca_c, 0)
#    Claw2, dummy, Cmaw2, xcp = lookup_data(0.087, ca_c, 0)
#    Clah1, Cdah0, Cmah1, xcp = lookup_data(Q_("0 deg"), ce_c, 0)
#    Clah2, dummy, Cmah2, xcp = lookup_data(0.087, ce_c, 0)
#    Cldh1, dummy, Cmdh1, xcp = lookup_data(0, ce_c, Q_("0 deg"))
#    Cldh2, dummy, Cmdh2, xcp = lookup_data(0, ce_c, 10)
#    
#    print(lin_ran_de)
#    Cl_alpha_w = (Claw2 - Claw1)/alpha_rad
#    Cl_alpha_h = (Clah2-Clah1)/alpha_rad
#    Cl_de_h = (Cldh2-Cldh1)/lin_ran_de
#    print(Cl_de_h)
#    Cm_alpha_w = (Cmaw2 - Cmaw1)/alpha_rad
#    Cm_alpha_h = (Cmah2 - Cmah1)/alpha_rad
#    Cm_de_h = (Cmdh2 - Cmdh1)/lin_ran_de
#    
#    q_dp = 0.5*rho*V**2                                       #Dynamic pressure
#    trim_mat1 = np.matrix([[(((Cdaw0- Cl_alpha_w)*S_w+(Cdah0 - Cl_alpha_h)*S_h)*q_dp).magnitude,
#                            (Cl_de_h * S_h * q_dp).magnitude],
#                            [(Cl_alpha_w * q_dp * S_w * (x_w-xcg) + Cl_alpha_h * q_dp *
#                             S_h * (x_h - xcg) + Cm_alpha_w * q_dp * S_w * MAC +
#                             Cm_alpha_h * q_dp * S_h * MAC_htail).magnitude,
#                             (Cl_de_h * q_dp * S_h * (x_h - xcg) +
#                             Cm_de_h * q_dp * S_h * MAC_htail).magnitude]])
#    trim_mat2 = np.matrix([[-W.magnitude],
#                           [0]])
#    trim_cond = np.linalg.solve(trim_mat1, trim_mat2)
#    alpha_trim = trim_cond[0,0] * Q_("rad")
#    de_trim = trim_cond[1,0] * Q_("deg")
#    Cltrimw, Cdtrimw, Cmtrimw, xcp = lookup_data(alpha_trim, ca_c, 0)
#    Cltrimh, Cdtrimh, Cmtrimh, xcp = lookup_data(alpha_trim, ce_c, de_trim)
#    Treq = (np.cos(alpha_trim)*(Cdtrimw*S_w + Cdtrimh * S_h)+np.sin(alpha_trim)*
#            (Cltrimw * S_w + Cltrimh * S_h)) * q_dp
#    print((((Cdaw0- Cl_alpha_w)*S_w+(Cdah0 - Cl_alpha_h)*S_h)*q_dp)*trim_cond[0] + 
#                            (Cl_de_h * S_h * q_dp)* (de_trim))
#    print(-W)
#    return alpha_trim, de_trim, Treq
# For roll
#de = -5.9
#da = 29.5

alpha_w_l = alpha_nose
alpha_w_r = alpha_nose
alpha_h = alpha_nose
V_local = V_inf

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

"============================================================================="
# ROLL
#for t_current in tlst:
#    #Wing
#    for i in range(0, len(kwlst)-1):
#        da_local = disc_wing_w[i,0]                        # Local aileron deflection of piece
#        b1 = kwlst[i]                                       # Y boundary left of piece
#        b2 = kwlst[i+1]                                     # Y boundary right of piece
#        y_i = (b1 + b2) / 2                                 # Y centre of piece
#        c1 = local_chord(abs(b1), c_r_w, c_t_w, half_b_w)   # Chord left of piece
#        c2 = local_chord(abs(b2), c_r_w, c_t_w, half_b_w)   # Chord right of piece
#        cc = local_chord(abs(y_i), c_r_w, c_t_w, half_b_w)  # Chord centre of piece
#        ca_c = (l_a / cc).magnitude                         # percentage aileron chord over local aileron
#        Sloc = (c1 + c2) / 2 * (b2 - b1)                    # Surface area of piece
#        
#        roll_induced_alpha = p * y_i / V_local
#        alpha_w = alpha_nose + roll_induced_alpha 
#        
#        Cl, Cd, Cm, xcp = lookup_data(alpha_w, ca_c, da_local)
#        Cn_w = -Cl*m.cos(alpha_w) - Cd*m.sin(alpha_w)
#        Ct_w = Cl*m.sin(alpha_w) - Cd*m.cos(alpha_w) 
#        
#        Fn_w_d = 0.5 * rho * V_local ** 2 * Sloc * Cn_w
#        Ft_w_d = 0.5 * rho * V_local ** 2 * Sloc * Ct_w
#        
#        disc_wing_w[i,1] = Fn_w_d.magnitude
#        disc_wing_w[i,2] = Ft_w_d.magnitude
#        disc_wing_w[i,3] = (Fn_w_d*y_i).magnitude
#        
#
#    # Wing simple
#    
#    dx_w = (y_mac+0.25*MAC)-xcg
#    Cl, Cd, Cm, xcp = lookup_data(alpha_nose, 0.2,0.)
#    Cn_w = -Cl*m.cos(alpha_nose) - Cd*m.sin(alpha_nose)
#    Ct_w = Cl*m.sin(alpha_nose) - Cd*m.cos(alpha_nose) 
#    
#    Fn_w = 0.5 * rho * V_local ** 2 * S_w * Cn_w
#    Ft_w = 0.5 * rho * V_local ** 2 * S_w * Ct_w
#
#    #Horizontal Tail
#    dx_h = l_h
#    Cl, Cd, Cm, xcp = lookup_data(alpha_h, 0.5,de)
#    Cn_h = -Cl*m.cos(alpha_h) - Cd*m.sin(alpha_h)
#    Ct_h = Cl*m.sin(alpha_h) - Cd*m.cos(alpha_h) 
#    
#    Fn_h = 0.5 * rho * V_local ** 2 * S_h * Cn_h
#    Ft_h = 0.5 * rho * V_local ** 2 * S_h * Ct_h
#    
#    # Moments and Forces
#    sum_Fn_w = sum(disc_wing_w[:,1])*Q_("N")
#    My = dx_w * sum_Fn_w + dx_h * Fn_h
#    Mx = sum(disc_wing_w[:,3])*Q_("N*m")
#    
#    p_dot = Mx/I_xx
#    q_dot = My/I_yy
#    
#    p += p_dot * dt
#    q += q_dot * dt
#    
#    plst.append(p.to('degree/s'))
#
#plt.plot(tlst,plst)
#print("Fz:",Fn_w+Fn_h+W)
"============================================================================="
# Wing
#de = -2.55
#Thrust = Q_("530 N")
#theta = alpha_nose
#q = Q_("0. rad/s")
#u = m.cos(alpha_nose) * V_inf
#w = m.sin(alpha_nose) * V_inf
#
#dx_w = (y_mac+0.25*MAC)-xcg 
#dx_h = l_h-xcg + Q_('0 m')
#
#qlst = []
#thetalst = []
#alst = []
#tlst = np.arange(0,5,dt.magnitude)
#Vlst = []
#
#for t_current in tlst:
#    if t_current >2:
#        de = -25
#    #Main Wing
#    alpha_w = alpha_nose + q*dx_w/V_local
#    Cl, Cd, Cm, xcp = lookup_data(alpha_w, 0.2,0.)
#    Cn_w = -Cl*m.cos(alpha_w) - Cd*m.sin(alpha_w)
#    Ct_w = Cl*m.sin(alpha_w) - Cd*m.cos(alpha_w) 
#    Fn_w = 0.5 * rho * V_local ** 2 * S_w * Cn_w
#    Ft_w = 0.5 * rho * V_local ** 2 * S_w * Ct_w
#    
#    #Horizontal Tail
#    alpha_h = alpha_nose + q*dx_h/V_local
#    Cl, Cd, Cm, xcp = lookup_data(alpha_h, ce_c,de)
#    Cn_h = -Cl*m.cos(alpha_h) - Cd*m.sin(alpha_h)
#    Ct_h = Cl*m.sin(alpha_h) - Cd*m.cos(alpha_h) 
#    
#    Fn_h = 0.5 * rho * V_local ** 2 * S_h * Cn_h
#    Ft_h = 0.5 * rho * V_local ** 2 * S_h * Ct_h
#    
#    Fx = Thrust + Ft_w + Ft_h - W * m.sin(theta)
#    Fz = Fn_w + Fn_h + W * m.cos(theta)
#    My = dx_w * Fn_w + dx_h * Fn_h
#    
#    u_dot = Fx/mtow - q*w
#    w_dot = Fz/mtow + q*u
#    q_dot = My/I_yy
#    
#    u += u_dot * dt
#    w += w_dot * dt
#    q += q_dot * dt
#    
#    theta += q * dt
#    alpha_nose = np.arctan(w/u) 
#
#    
#    V_local = np.sqrt(u**2+w**2)
#    
#    alst.append(alpha_nose.to('degree').magnitude)
#    qlst.append(q.to('degree /s').magnitude)
#    thetalst.append(theta.to('degree').magnitude)
#    
#    Vlst.append(V_local.magnitude)
#print('finished')
#
#plt.figure()
#plt.plot(tlst,qlst,label='q')
#plt.plot(tlst,thetalst,label='theta')
#plt.plot(tlst,alst,label='alpha')
#plt.legend()
##plt.figure()
##plt.plot(tlst,Vlst,label='velocity')
##plt.legend()
#plt.show()
"============================================================================="
#Vertical Tail
dr = 0
dx_v = l_v-xcg + Q_('1 m')
u = m.cos(beta_nose) * V_inf
v = m.sin(beta_nose) * V_inf
r = Q_('0. rad/s')
Thrust = Q_('1000 N')
alpha_v = alpha_nose
betalst = []
psilst = []
rlst = []
tlst = []
for t_current in np.arange(0,4,dt.magnitude):
    if t_current>1:
        dr = 25
        
    beta_v  = beta_nose - r*dx_v/V_inf
    
    Cl, Cd, Cm, xcp = lookup_data(-beta_v, cr_c,dr)
    Cn_v = -Cl*m.cos(beta_v) - Cd*m.sin(beta_v)
    Ct_v = Cl*m.sin(beta_v) - Cd*m.cos(beta_v) 
    
    Fn_v = 0.5 * rho * V_local ** 2 * S_v * Cn_v
    Ft_v = 0.5 * rho * V_local ** 2 * S_v * Ct_v
    
    Cl_w,Cd_w,dummy,dummy = lookup_data(alpha_nose,0.1,0)
    Cd_w = -Cl_w*m.cos(alpha_nose) - Cd_w*m.sin(alpha_nose)
    Ft_wh = 0.5 * rho * V_local**2 * Cd_w * (S_w + S_h)
    
    
    Fx = Thrust + Ft_wh + Ft_v  - W * m.sin(alpha_nose)
    Fz = Fn_v
    
    u_dot = Fx/mtow + r*v
    v_dot = Fz/mtow - r*u
    
    u += u_dot * dt
    v += v_dot * dt
    
    Mz = dx_v * Fn_v
    r_dot = Mz/I_zz
    
    r += r_dot * dt
    Psi += r * dt
    
    beta_nose = np.arctan(v/u)
    
    betalst.append(m.degrees(beta_nose))
    psilst.append(m.degrees(Psi))
    rlst.append(r.to('degree/s').magnitude)
    tlst.append(t_current)
plt.plot(tlst, betalst, label='beta')
plt.plot(tlst, psilst, label='psi')
plt.plot(tlst,rlst,label='r')
plt.legend()
plt.show()

