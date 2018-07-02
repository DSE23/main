"""
Name: Gyroscopic simulation
Department: Propulsion and Aircraft Systems
Last updated: 02/07/2018 15:06 by Ties
"""

"""
This file simulates gyroscopic time responses
"""


from APMonitor.apm import *
import numpy as np
import matplotlib.pyplot as plt

n = 501
tiime = 5
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

# fig, ax = plt.subplots()
# ax.plot(t, phi_dot, label="Phi_dot")
# ax.plot(t, theta_dot, label="Theta_dot")
# ax.plot(t, psi_dot, label="Psi_dot")
# plt.ylabel('Angles [rad/s]')
# plt.xlabel('Time [s]')
# plt.title('Angles with 1500Nm elevator moment')
# plt.legend()
# plt.show()

# print(solution1)
fig, ax = plt.subplots()
ax.plot(t, p, 'k', label="p")
ax.plot(t, q, 'k--', label="q")
ax.plot(t, r, 'k-.', label="r")
plt.ylabel('Rate [rad/s]')
plt.xlabel('Time [s]')
# plt.title('Response with 5000Nm rudder moment')
plt.legend()
plt.show()

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