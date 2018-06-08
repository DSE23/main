"""
Name: Engine Mounts
Department: Propulsion and Aircraft Systems
Last updated: 05/06/2018 16:52 by Ties
"""

"""
This file calculates the maximum forces experienced by the engine mounts, 
for use by the structures department
"""

import sys
import scipy as sp
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

# Import data
import Geometry
import Engine
import Propeller
import Firewall

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, pointing towards tail
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_v_dist(inp):
    global v_dist
    v_dist = inp


def initialise_h_dist(inp):
    global h_dist
    h_dist = inp


def initialise_length(inp):
    global length
    length = inp


def initialise_mass(inp):
    global mass
    mass = inp


def initialise_xcg(inp):
    global xcg
    xcg = inp


# End defining global variables

# Start assigning values to variables
# Vertical distance between engine mount attachment points to firewall
v_dist = Q_("1 m")  # DUMMY, can be chosen
v_dist.ito(ureg.m)
# Horizontal distance between engine mount attachment points to firewall
h_dist = Q_("1 m")  # DUMMY, can be chosen
h_dist.ito(ureg.m)
# Distance between engine mounts and firewall
length = Firewall.eng_clearance  # From Firewall file
length.ito(ureg.m)
# Engine mounts mass
mass = Q_("5 kg")  # DUMMY
# Engine mounts cg, assumed at end of engine (engine attachments are not at the end of the engine)
xcg = Engine.length

# Start calculations

# Load case
# Max load factor (defined from requirements):
loadfactor = 20
# Gravitational acceleration
gravity = Q_("9.81 m/s**2")
# Roll acceleration, assumed about crankshaft
roll_acc = Q_("1 rad/s**2")  # DUMMY

# Determine max engine vertical force (in Newtons)
f_y_eng = Engine.mass * gravity * loadfactor
f_y_eng.ito(ureg.newton)

# Determine max prop vertical force (in Newtons)
f_y_prop = Propeller.mass * gravity * loadfactor
f_y_prop.ito(ureg.newton)

# Determine max mount vertical force (in Newtons)
f_y_mount = mass * gravity * loadfactor
f_y_mount.ito(ureg.newton)

# Determine reaction forces
# Sum of forces in y-direction, r_y is engine mount reaction force
r_y_total = f_y_prop+f_y_eng

# Moments
# Component moments about axis parallel to Z-axis through bottom engine mount
m_z_eng = f_y_eng * (Firewall.xcg - Engine.xcg)
m_z_prop = f_y_prop * (Firewall.xcg - Propeller.xcg)
m_z_mount = f_y_mount * (Firewall.xcg - xcg)
# Sum of moments about axis parallel to Z-axis through bottom engine mount
r_x_1 = (m_z_eng + m_z_prop) / v_dist
# Sum of forces in x-direction
r_x_2 = -r_x_1

print("Vertical shear force per mount: {}".format(r_y_total/4))
print("Normal force per top mount: {}".format(r_x_1/2))
print("Normal force per bottom mount: {}".format(r_x_2/2))
