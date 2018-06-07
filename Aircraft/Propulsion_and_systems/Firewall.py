"""
Name: Firewall
Department: Propulsion and Aircraft Systems
Last updated: 07/06/2018 12:55 by Ties
"""

import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

import Engine

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, pointing towards tail
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_eng_clearance(inp):
    global eng_clearance
    eng_clearance = inp


def initialise_xcg(inp):
    global xcg
    xcg = inp


# End defining global variables

# Start assigning values to variables
# Distance between engine and firewall, assumed to be 10cm
eng_clearance = Q_("10 cm")
# Position behind datum
xcg = Engine.length + eng_clearance
