"""
Name: Engine Mounts
Department: Propulsion and Aircraft Systems
Last updated: 05/06/2018 16:52 by Ties
"""

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry from the Misc folder

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

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, backwards
# Y axis: up
# Z axis: left

# Engine dry mass as provided by Lycoming
enginedrymass = Q_('446 lbs')
# Assume factor of 10% to achieve wet mass
enginewetmass = 1.1*enginedrymass

# Engine mass moments of inertia about engine cg
engine_ixg = Q_("84.4 inch*lbf*s**2")
engine_iyg = Q_('93.5 inch*lbf*s**2')
engine_izg = Q_('145.8 inch*lbf*s**2')

# Engine dimensions
enginelength = Q_('39.34 in')
enginewidth = Q_('34.25 in')
engineheight = Q_('26.46 in')

# Distance from engine cg to origin
