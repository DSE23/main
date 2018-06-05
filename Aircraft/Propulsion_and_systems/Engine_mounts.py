"""
Name: Engine Mounts
Department: Propulsion and Aircraft Systems
Last updated: 05/06/2018 16:52 by Ties
"""

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder


def initialise_enginemass(inp):
    global enginedrymass
    enginedrymass = inp


def initialise_engineIxg(inp):
    global engineIxg
    engineIxg = inp


def initialise_engineIzg(inp):
    global engineIzg
    engineIzg = inp


def initialise_engineIyg(inp):
    global engineIyg
    engineIyg = inp


def initialise_enginelength(inp):
    global enginelength
    enginelength = inp


def initialise_enginewidth(inp):
    global enginewidth
    enginewidth = inp


def initialise_engineheight(inp):
    global engineheight
    engineheight = inp


# Engine dry mass as provided by Lycoming
enginedrymass = Q_('446 lbs')
# Assume factor of 10% to achieve wet mass
enginewetmass = 1.1*enginedrymass

# Mass moments of inertia
# About the axis parallel to the crankshaft centerline
engineIxg = Q_("84.4 inch*lbf*s**2")
# About the vertical axis
engineIyg = Q_('93.5 inch*lbf*s**2')
# About the axis parallel to the centerline
engineIzg = Q_('145.8 inch*lbf*s**2')

#Engine dimensions
enginelength = Q_('39.34 in')
enginewidth = Q_('34.25 in')
engineheight = Q_('26.46 in')

