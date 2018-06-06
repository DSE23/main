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