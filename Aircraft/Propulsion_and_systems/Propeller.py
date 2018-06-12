"""
Name: Propeller
Department: Propulsion and Aircraft Systems
Last updated: 06/06/2018 10:25 by Ties
"""

import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, pointing towards tail
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_mass(inp):
    global mass
    mass = inp


def initialise_propxcg(inp):
    global xcg
    xcg = inp


def initialise_propinertia(inp):
    global ixg
    ixg = inp


# End defining global variables

# Start assigning values to variables
# Propeller mass
mass = Q_("30 kg")  # Based on MT-propeller (4-bladed)

# CG
xcg = Q_("-20 cm")  # DUMMY VALUE, NOT KNOWN YET, negative because in front of datum

# Mass moment of inertia
ixg = Q_("0.9 kg * m**2")  # Based on slightly different propeller
