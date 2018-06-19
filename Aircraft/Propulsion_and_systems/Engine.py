"""
Name: Engine
Department: Propulsion and Aircraft Systems
Last updated: 06/06/2018 10:21 by Ties
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


def initialise_mass(inp):
    global mass
    mass = inp


def initialise_ixg(inp):
    global ixg
    ixg = inp


def initialise_izg(inp):
    global izg
    izg = inp


def initialise_iyg(inp):
    global iyg
    iyg = inp


def initialise_length(inp):
    global length
    length = inp


def initialise_width(inp):
    global width
    width = inp


def initialise_height(inp):
    global height
    height = inp


def initialise_xcg(inp):
    global xcg
    xcg = inp


def initialise_ycg(inp):
    global ycg
    ycg = inp


def initialise_zcg(inp):
    global zcg
    zcg = inp

# End defining global variables

# Start assigning values to variables
# Engine dry mass as provided by Lycoming
drymass = Q_("446 lbs")
drymass.ito(ureg.kg)
# Assume additional 10% to include fluids and hoses
mass = 1.1*drymass

# Engine mass moments of inertia about engine cg, as provided by Lycoming
ixg = Q_("84.4 inch*lbf*s**2")
iyg = Q_("93.5 inch*lbf*s**2")
izg = Q_("145.8 inch*lbf*s**2")
# Transfer to SI units
ixg.ito(ureg("kg*m**2"))
iyg.ito(ureg("kg*m**2"))
izg.ito(ureg("kg*m**2"))

# Engine dimensions
length = Q_("39.34 in")
width = Q_("34.25 in")
height = Q_("26.46 in")
# Transfer to SI units
length.ito(ureg.m)
width.ito(ureg.m)
height.ito(ureg.m)

# Engine cg location
# Assumed cg is on crankshaft & in the middle of the engine length
xcg = length/2
ycg = Q_("0 m")
zcg = Q_("0 m")
