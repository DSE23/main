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

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, backwards
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


# Test

# End defining global variables

# Start assigning values to variables
# Vertical distance between engine mount attachment points to firewall
v_dist = Q_("10 m")  # DUMMY
v_dist.ito(ureg.m)
# Horizontal distance between engine mount attachment points to firewall
h_dist = Q_("10 m")  # DUMMY
h_dist.ito(ureg.m)
# Distance between engine mounts and firewall
length = Q_("10 m")  # DUMMY
length.ito(ureg.m)

# Load case
# Max load factor (defined from requirements):
loadfactor = 20
# Gravitational acceleration
gravity = Q_("9.81 m/s**2")

# Determine max engine vertical force (in Newtons)
f_y_eng = Engine.mass * gravity * loadfactor
f_y_eng.ito(ureg.newton)

# Determine max prop vertical force (in Newtons)
f_y_prop = Propeller.mass * gravity * loadfactor
f_y_prop.ito(ureg.newton)

# Calculate x-distance from engine cg to firewall
d_eng_firewall = firewall_xcg - Engine.xcg

# Determine reaction forces
# Sum of forces in y-direction, R_y is engine mount reaction force
r_y_total = f_y_prop+f_y_eng

# Sum of moments about axis parallel to Z-axis through bottom engine mount
r_x_1 = (f_y_eng * (Engine.length-Engine.xcg) + f_y_prop * (Engine.length-Propeller.xcg)) / v_dist
print(r_x_1)
