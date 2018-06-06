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
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

import Geometry

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, backwards
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_propmass(inp):
    global propmass
    propmass = inp


def initialise_propcg(inp):
    global propcg
    propcg = inp


# End defining global variables

# Propeller mass
propmass = Q_("10000 kg") #DUMMY VALUE, NOT KNOWN YET

# Propeller CG
propcg = Q_("-1000 cm") #DUMMY VALUE, NOT KNOWN YET

# Import aircraft cg location from geometry
