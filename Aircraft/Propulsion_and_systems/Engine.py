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


def initialise_enginemass(inp):
    global enginedrymass
    enginedrymass = inp


def initialise_engine_ixg(inp):
    global engine_ixg
    engine_ixg = inp


def initialise_engine_izg(inp):
    global engine_izg
    engine_izg = inp


def initialise_engine_iyg(inp):
    global engine_iyg
    engine_iyg = inp


def initialise_enginelength(inp):
    global enginelength
    enginelength = inp


def initialise_enginewidth(inp):
    global enginewidth
    enginewidth = inp


def initialise_engineheight(inp):
    global engineheight
    engineheight = inp


def initialise_engine_xcg(inp):
    global engine_xcg
    engine_xcg = inp


def initialise_engine_ycg(inp):
    global engine_ycg
    engine_ycg = inp


def initialise_engine_zcg(inp):
    global engine_zcg
    engine_zcg = inp


# End defining global variables

# Start assigning values to variables
# Engine dry mass as provided by Lycoming
enginedrymass = Q_("446 lbs")
# Assume factor of 10% to achieve wet mass
enginewetmass = 1.1*enginedrymass

# Engine mass moments of inertia about engine cg
engine_ixg = Q_("84.4 inch*lbf*s**2")
engine_iyg = Q_("93.5 inch*lbf*s**2")
engine_izg = Q_("145.8 inch*lbf*s**2")

# Engine dimensions
enginelength = Q_("39.34 in")
enginewidth = Q_("34.25 in")
engineheight = Q_("26.46 in")

# Engine cg location
# Assumed cg is on crankshaft & in the middle of the engine length
engine_xcg = enginelength/2
engine_ycg = 0
engine_zcg = 0
