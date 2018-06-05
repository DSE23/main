"""
Name: Engine Mounts
Department: Propulsion and Aircraft Systems
Last updated: 05/06/2018 16:52 by Ties
"""

from pint import UnitRegistry

unit = UnitRegistry()


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
enginedrymass = unit('446 lbs')
# Assume factor of 10% to achieve wet mass
enginewetmass = 1.1*enginedrymass

# Mass moments of inertia
# About the axis parallel to the crankshaft centerline
engineIxg = unit("84.4 inch*lbf*s**2")
# About the vertical axis
engineIyg = unit('93.5 inch*lbf*s**2')
# About the axis parallel to the centerline
engineIzg = unit('145.8 inch*lbf*s**2')

#Engine dimensions
enginelength = unit('39.34 in')
enginewidth = unit('34.25 in')
engineheight = unit('26.46 in')

