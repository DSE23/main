"""
Name: Gyroscopic effects
Department: Propulsion and Aircraft Systems
Last updated: 15/06/2018 13:58 by Ties
"""

"""
This file calculates gyroscopic effects of the propeller during turns and when applying a moment
"""

import sys
import scipy as sp
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

from Propulsion_and_systems import Propeller
from APMonitor.apm import *
import numpy as np
import matplotlib.pyplot as plt

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centreline, pointing towards tail
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_rpm(inp):
    global rpm
    rpm = inp


# Start assigning values
rpm = Q_("2700 rpm")
rpm.ito(ureg("rad/min"))

# Make function which calculates gyro accelerations (prop only!) from moment and rate inputs


def input_moment(yaw_moment, pitch_moment, yaw_rate, pitch_rate):
    # Make compatible for inputs with and without units
    if yaw_moment.__class__ == int or yaw_moment.__class__ == float:
        yaw_moment *= ureg("newton*meter")
    if pitch_moment.__class__ == int or pitch_moment.__class__ == float:
        pitch_moment *= ureg("newton*meter")
    if yaw_rate.__class__ == int or yaw_rate.__class__ == float:
        yaw_rate *= ureg("rad/s")
    if pitch_rate.__class__ == int or pitch_rate.__class__ == float:
        pitch_rate *= ureg("rad/s")

    omega_y_dot = (yaw_moment - (Propeller.ixg - Propeller.izg) * pitch_rate * rpm) / Propeller.izg
    omega_z_dot = (pitch_moment - (Propeller.izg - Propeller.ixg) * yaw_rate * rpm) / Propeller.izg

    omega_y_dot.ito(ureg("rad/s**2"))
    omega_z_dot.ito(ureg("rad/s**2"))

    print("Yaw acceleration: {}".format(omega_y_dot))
    print("Pitch acceleration: {}".format(omega_z_dot))
    return [omega_y_dot, omega_z_dot]


# Test the function!
# input_moment(0, Q_("1 N*m"), 0, 0)

# Make function which calculates gyro moments from rate and acceleration inputs


def input_acceleration(yaw_acceleration, pitch_acceleration, yaw_rate, pitch_rate):
    # Make compatible for inputs with and without units
    if yaw_acceleration .__class__ == int or yaw_acceleration.__class__ == float:
        yaw_acceleration *= ureg("rad/s**2")
    if pitch_acceleration .__class__ == int or pitch_acceleration.__class__ == float:
        pitch_acceleration *= ureg("rad/s**2")
    if yaw_rate.__class__ == int or yaw_rate.__class__ == float:
        yaw_rate *= ureg("rad/s")
    if pitch_rate.__class__ == int or pitch_rate.__class__ == float:
        pitch_rate *= ureg("rad/s")

    m_y = Propeller.izg * yaw_acceleration + (Propeller.ixg-Propeller.izg) * pitch_rate * rpm
    m_z = Propeller.izg * pitch_acceleration + (Propeller.izg - Propeller.ixg) * yaw_rate * rpm

    m_y.ito(ureg("newton * meter"))
    m_z.ito(ureg("newton * meter"))

    print("Yaw moment: {}".format(m_y))
    print("Pitch moment: {}".format(m_z))
    return [m_y, m_z]


# Test the function!
input_acceleration(0, 0, 0, Q_("60 deg/s"))

y_acc = input_moment(0, Q_("1 N*m"), 0, 0)[0]
z_acc = input_moment(0, Q_("1 N*m"), 0, 0)[1]
input_acceleration(y_acc, z_acc, 0, 0)


# Solver shit

n = 501
tiime = 15
time = np.linspace(0, tiime, n)
p = np.ones(n)
x = 2 * np.ones(n)

fid = open('data.csv', 'w')
# write time row
fid.write('time, ')
for i in range(n-1):
    fid.write(str(time[i]) +  ', ')
fid.write(str(time[n-1]) + '\n')

# write 'p' row (input parameter)
fid.write('p, ')
for i in range(n-1):
    fid.write(str(p[i]) +  ', ')
fid.write(str(p[n-1]) + '\n')

# write 'x' row (state variable initialization)
fid.write('x, ')
# imode: http://apmonitor.com/wiki/index.php/Main/Modes
# for imode=4-6, include all initialization values
# for imode=7-9, include only the initial condition for variables
imode = 7
if ((imode>=4) and (imode<=6)):
    for i in range(n-1):
        fid.write(str(x[i]) + ', ')
    fid.write(str(x[n-1]) + '\n')
else:
    fid.write(str(x[0]) + ', ')
    for i in range(1,n-1):
        fid.write('-, ')
    fid.write('-\n')

# close file
fid.close()

s = 'http://byu.apmonitor.com'
a = 'gyro'

apm(s, a, 'clear all')
apm_load(s, a, 'gyro.apm')
csv_load(s, a, 'data.csv')

apm_option(s, a, 'apm.time_shift', 1)
apm_option(s, a, 'apm.imode', imode)
output1 = apm(s, a, 'solve')

solution1 = apm_sol(s, a)

t = np.array(solution1['time'])
dt = t[1]
p = -np.array(solution1['gyro.omega_x_a'])
q = -np.array(solution1['gyro.omega_z'])
r = -np.array(solution1['gyro.omega_y'])

phi_dot = np.zeros(n)
theta_dot = np.zeros(n)
psi_dot = np.zeros(n)

phi = np.zeros(n)
theta = np.zeros(n)
psi = np.zeros(n)

for i in range(n-1):
    phi_dot[i] = p[i] + np.sin(phi[i])*np.tan(theta[i]) * q[i] + np.cos(phi[i])*np.tan(theta[i])*r[i]
    theta_dot[i] = np.cos(phi[i]) * q[i] - np.sin(phi[i]) * r[i]
    psi_dot[i] = (np.sin(phi[i])/np.cos(theta[i]))*q[i] + (np.cos(phi[i])/np.cos(theta[i]))*r[i]

    phi[i+1] = phi[i] + phi_dot[i] * dt
    theta[i+1] = theta[i] + theta_dot[i] * dt
    psi[i+1] = psi[i] + psi_dot[i] * dt

fig, ax = plt.subplots()
ax.plot(t, phi, label="Phi")
ax.plot(t, theta, label="Theta")
ax.plot(t, psi, label="Psi")
plt.ylabel('Angles [rad/s]')
plt.xlabel('Time [s]')
plt.title('Angles with 5000Nm rudder moment')
plt.legend()
plt.show()

# print(solution1)
# fig, ax = plt.subplots()
# ax.plot(solution1['time'], solution1['gyro.omega_y'], label="Yaw")
# ax.plot(solution1['time'], solution1['gyro.omega_z'], label="Pitch")
# ax.plot(solution1['time'], solution1['gyro.omega_x_a'], label="Roll")
# plt.ylabel('Rates [rad/s]')
# plt.xlabel('Time [s]')
# plt.title('Rates with 5000Nm rudder moment')
# plt.legend()
# plt.show()

# apm.coldstart = 1
# apm.ctrl_units = 1
# apm.ctrl_time = 0.01
# apm.pred_time = 1000
# apm.pred_hor = 1000
# apm.hist_units = 1
# apm.hist_hor = 1000
# apm.web = 1
# y = apm_solve("gyro", 7, "http://byu.apmonitor.com")
# #z = apm_web("http://byu.apmonitor.com", 'gyro')
#
# print(y)
#omega_y = y['omega_y']
#omega_z = y['omega_z']
