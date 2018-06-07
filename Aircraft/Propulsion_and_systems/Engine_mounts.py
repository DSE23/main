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

# Test

# End defining global variables

# Start assigning values to variables
# Vertical distance between engine mounts
v_dist = Q_("20 in")  # DUMMY
v_dist.ito(ureg.m)
# Horizontal distance between engine mounts
h_dist = Q_("10000000 m")  # DUMMY
h_dist.ito(ureg.m)

# Load case
# Max load factor (defined from requirements):
loadfactor = 20
# Gravitational acceleration
gravity = Q_("9.81 m/s**2")

# Determine max engine vertical force (in Newtons)
F_y_eng = Engine.mass * gravity * loadfactor
F_y_eng.ito(ureg.newton)

# Determine max prop vertical force (in Newtons)
F_y_prop = Propeller.mass * gravity * loadfactor
F_y_prop.ito(ureg.newton)

# Determine moment about Z axis at mounts
# M_z_eng = F_y_eng * (Engine.length-Engine.xcg)
# M_z_prop = F_y_prop * (Engine.length-Propeller.xcg)

# Determine reaction forces
# Sum of forces in y-direction, R_y is engine mount reaction force
R_y_total = F_y_prop+F_y_eng

# Sum of moments about axis parallel to Z-axis through bottom engine mount
R_x_1 = (F_y_eng * (Engine.length-Engine.xcg) + F_y_prop * (Engine.length-Propeller.xcg)) / v_dist
print(R_x_1)
