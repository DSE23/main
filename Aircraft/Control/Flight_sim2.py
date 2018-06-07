"""
Created on Wed May 18 12:15:59 2016

@author: Jordy
"""
import math
import time
import Calc_ISA_100km
import cl_cd

#Init sim
labda = 51.95203                # [deg] latitude
mu    = 4.43293                 # [deg] longitude

theta = 0.                      # [rad] angle around y-axis
phi   = 0.                      # [rad] angle around x-axis
psy   = 0.                      # [rad] angle around z-axis
M = 0.                          # Moment around y-axis
L = 0.                          # Moment around x-axis
N = 0.                          # Moment around z-axis

q = 0.                          # [rad/s] angular velocity around y-axis
q_dot = 0.                      # [rad/s2] angular acceleration around y-axis
p = 0.                          # [rad/s] angular velocity around x-axis
p_dot = 0.                      # [rad/s2] angular acceleration around x-axis
r = 0.                          # [rad/s] angular velocity around z-axis
r_dot = 0.                      # [rad/s2] angular acceleration around z-axis

v = 0.                          # [m/s] linear velocity along y-axis
v_dot = 0.                      # [m/s2] linear acceleration along y-axis
u = 90.                         # [m/s] linear velocity along x-axis 
u_dot = 0.                      # [m/s2] linear acceleration along x-axis
w = 0.                          # [m/s] linear velocity along z-axis
w_dot = 0.                      # [m/s2] linear acceleration along z-axis

X = 0.
Z = 0.
Y = 0.

h = 500.                        # [m] altitude

e0  = 1.                        # [-] Quaternion parameter
e1  = 0.                        # [-] Quaternion parameter
e2  = 0.                        # [-] Quaternion parameter
e3  = 0.                        # [-] Quaternion parameter

a11 = e0*e0 + e1*e1 - e2*e2 - e3*e3     # [-] DCM parameter
a12 = 2 * (e1*e2 - e0*e3)               # [-] DCM parameter
a13 = 2 * (e0*e2 + e1*e3)               # [-] DCM parameter
a21 = 2 * (e1*e2 + e0*e3)               # [-] DCM parameter
a22 = e0*e0 - e1*e1 + e2*e2 - e3*e3     # [-] DCM parameter
a23 = 2 * (e2*e3 - e0*e1)               # [-] DCM parameter
a31 = 2 * (e1*e3 - e0*e2)               # [-] DCM parameter
a32 = 2 * (e2*e3 + e0*e1)               # [-] DCM parameter
a33 = e0*e0 - e1*e1 - e2*e2 + e3*e3     # [-] DCM parameter

# Constants
Rearth= 6371000                 # [m] Earth radius
g     = 9.80665                 # [m/s2] gravitational acceleration
demax = 30.                     # [deg] maximum elevator deflection (downwwards)
demin = -30.                    # [deg] minimum elevator deflection (upwards)

# Here comes StefX into play
Ixx   = 1089.319031             # [kg*m2] roll inertia
Iyy   = 1564.932247             # [kg*m2] pitch inertia
Izz   = 2629.662716             # [kg*m2] yaw inertia
Ixz   = 0.                      # [kg*m2] lateral cross inertia
                 
m     = 889.3844425             # [kg] aircraft mass
W     = m*g                     # [N] aircraft weight
S     = 11.74274748             # [m2] wing surface area
c     = 1.53                    # [m] wing mean aerodynamic chord
b     = 8.04                    # [m] wing span
cg    = 0.14                    # [m] Centre of gravity location (1.453765363m from nose)



Cybeta= -0.14                   # [1/rad]
Cydr  = 5.57                    # [1/rad]           REALLY HIGH
Clbeta= 71*10**(-5)             # [1/rad]
Clp   = -0.63                   # [1/rad]
Clr   = 0.0383                  # [1/rad]
Clda  = 3.92                    # [1/rad]
Cldr  = 0.0147##                  # [1/rad]
Cm0   = 0.                      # [1/rad]           CHECK!
Cma   = -0.88                   # [1/rad]
Cmadot= -3.14                   # [1/rad]
Cmq   = -4.28                   # [1/rad]
Cmde  = -1.28##                   # [1/rad]
Cnbeta= 0.07729                 # [1/rad]
Cnp   = 0.0028                  # [1/rad]           SHOULD BE NEGATIVE?
Cnr   = -0.08961                # [1/rad]
Cnda  = -0.053##                  # [1/rad]
Cndr  = -0.0657##                 # [1/rad]


#Init clock
t = time.clock()
t0 = t
tstart = t
Thrust = 0.
running = True
# Start Simulation
while running:
    t  = time.clock()
    dt = t-t0 
    t0 = t
    trunning = t-tstart
    
    # 1) Input by joystick
    de = 0.
    da = 0.
    dr = 0.
    dp = 0.
    df = 0.
    
    # 2) Compute the angles of attack
    Vc = math.sqrt(u*u + v*v + w*w)
    
    aoa = math.atan2(w,u)
    aoa_dot = (u*w_dot - w*u_dot) / math.sqrt(u*u + w*w)

    beta = (v_dot * math.sqrt(u*u + w*w) - v*(u*u_dot + w*w_dot))/(math.sqrt(u*u+w*w) * (u*u + v*v +w*w))
    
    # 3) Coefficients of aerodynamic forces
    Cl = cl_cd.nonlinear(math.degrees(aoa))[0]
    Cd = cl_cd.nonlinear(math.degrees(aoa))[1]
    
    # 4) Coefficients of aerodynamic moments in pitch
    # 5) Coefficients of aerodynamic moments in roll
    # 6) Coefficients of aerodynamic moments in yaw
    
    
    # 7) Body frame forces
    Lift = 0.5 * Calc_ISA_100km.isacal(h)[2] * Vc*Vc * Cl * S
    Drag = 0.5 * Calc_ISA_100km.isacal(h)[2] * Vc*Vc * Cd * S
    sideforce = 0.5 * Calc_ISA_100km.isacal(h)[2] * Vc*Vc * S * (Cydr*dr + Cybeta*beta)
    
    # 8) Engine forces and moments
    # 9) Gear forces and moments
    # 10) Resolve body frame forces
    Fx = Lift * math.sin(aoa) - Drag * math.cos(aoa) - m*g * math.sin(theta) + Thrust
    Fy = m*g*math.sin(phi)*math.cos(theta)
    Fz = - Lift * math.cos(aoa) - Drag * math.sin(aoa) + m*g * math.cos(theta) * math.cos(phi)
    
    # 11) Body frame accelerations
    u_dot = Fx/m
    v_dot = Fy/m
    w_dot = Fz/m
    
    # 12) Body frame aerodynamic velocities    
    u = u + u_dot * dt
    v = v + v_dot * dt
    w = w + w_dot * dt
    
    # 13) Wind components
    # 14) Turbulance componenents
    
    # 15) Earth velocities
    Vn = u*a11 + v*a12 + w*a13
    Ve = u*a21 + v*a22 + w*a23
    Vd = u*a31 + v*a32 + w*a33
    
    # 16) Latitude and Longitude rates
    labda_dot = Vn/(Rearth+h)
    mu_dot    = Ve/(math.cos(math.radians(labda)) * (Rearth+h))
    h_dot     = -Vd
    
    # 17) Aircraft position
    labda = labda + labda_dot * dt
    mu    = mu + mu_dot * dt
    h     = h  + h_dot * dt
    
    # 18) Body rates in stability axis
    p_stab = p * math.cos(aoa) + r * math.sin(aoa)
    r_stab = r * math.cos(aoa) - p * math.sin(aoa)
    
    # 19) Body frame moments in stability axis
    Mstab = 0.5 * Calc_ISA_100km.isacal(h)[2] * Vc*Vc * S * c*(Cm0 + Cma*aoa + Cmde*math.radians(de)) + 0.25 * Calc_ISA_100km.isacal(h)[2] * Vc * S * c*c * (Cmq*q + Cmadot * aoa_dot)
    Lstab = 0.5 * Calc_ISA_100km.isacal(h)[2] * Vc*Vc * S * b*(Clbeta*beta + Clda*math.radians(da) + Cldr*math.radians(dr)) + 0.25 * Calc_ISA_100km.isacal(h)[2] * Vc * S * b*b*(Clp*p_stab + Clr*r_stab)
    Rstab = 0.5 * Calc_ISA_100km.isacal(h)[2] * Vc*Vc * S * b*(Cnbeta*beta + Cnda*math.radians(da) + Cndr*math.radians(dr)) + 0.25 * Calc_ISA_100km.isacal(h)[2] * Vc * S * b*b*(Cnp*p_stab + Cnr*r_stab)
    
    # 20) Body frame moments in body frame axis
    M = Mstab + Lift*(cg-0.25)*c*math.cos(aoa) + Drag*(cg-0.25)*c*math.sin(aoa)
    L = Lstab*math.cos(aoa) - Rstab*math.cos(aoa) 
    R = Rstab*math.cos(aoa) + Lstab*math.sin(aoa) - sideforce*(cg-0.25)*c 
    
    # 21) Body frame angular accelerations
    p_dot = (L + (Iyy-Izz)*q*r + Ixz*(r_dot*p + p*q))/Ixx
    q_dot = (M + (Izz-Ixx)*r*p + Ixz*(r*r -p*p))/Iyy
    r_dot = (R + (Ixx-Iyy)*p*q + Ixz*(p_dot - q*r))/Izz
    
    # 22) Body rates
    p = p + p_dot * dt    
    q = q + q_dot * dt
    r = r + r_dot * dt
    
    # 23) quaternions
    labda_quat  = 1 - (e0*e0 + e1*e1 + e2*e2 + e3*e3)

    e0_dot = -0.5 * (e1*p + e2*q + e3*r) + labda_quat*e0
    e1_dot = 0.5 * (e0*p + e2*r - e3*q) + labda_quat*e1
    e2_dot = 0.5 * (e0*q + e3*p - e1*r) + labda_quat*e2
    e3_dot = 0.5 * (e0*r + e1*q - e2*p) + labda_quat*e3
    
    e0 = e0 + e0_dot * dt
    e1 = e1 + e1_dot * dt
    e2 = e2 + e2_dot * dt
    e3 = e3 + e3_dot * dt
    
    # 24) DCM
    a11 = e0*e0 + e1*e1 - (e2*e2) - (e3*e3)
    a12 = 2 * (e1*e2 - (e0*e3))
    a13 = 2 * (e0*e2 + e1*e3)
    a21 = 2 * (e1*e2 + e0*e3)
    a22 = e0*e0 - (e1*e1) + e2*e2 - (e3*e3)
    a23 = 2 * (e2*e3 - (e0*e1))
    a31 = 2 * (e1*e3 -( e0*e2))
    a32 = 2 * (e2*e3 + e0*e1)
    a33 = e0*e0 - (e1*e1) - (e2*e2) + e3*e3
    
    
    # 25) Euler angles
    theta = math.asin(-a31)
    phi   = math.atan2(a32,a33)
    psy   = math.atan2(a21,a11)
    
    
    if trunning > 1:
        running = False
        
    print ('Fx:',Fx)
    print ('Theta:',math.degrees(theta))
    print ("Lift Value:",Lift * math.sin(aoa))
    print ("Drag Value:",- Drag * math.cos(aoa))
    print ("Mass Value:",- m*g * math.sin(theta))
    print ("")
    Thrust = 2600.


print("Ready")
