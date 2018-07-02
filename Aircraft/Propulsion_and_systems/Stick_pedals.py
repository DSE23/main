"""
Name: Stick
Department: Propulsion and Aircraft Systems
Last updated: 12/06/2018 15:05 by Ties
"""

import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_stick_length(inp):
    global l_s_t
    l_s_t = inp


def initialise_bottom_length(inp):
    global l_s_b
    l_s_b = inp


def initialise_pedal_travel(inp):
    global d_p
    d_p = inp


def initialise_stick_defl_e(inp):
    global d_e
    d_e = inp


def initialise_stick_defl_r(inp):
    global d_r
    d_r = inp


# End defining global variables

# Start assigning values to variables
# Total stick length
l_s_t = Q_("81 cm")  # As defined by Gijs

# Stick length from hinge point to pilot's hand
l_s = Q_("60 cm")  # DUMMY

# Stick length below hinge point
l_s_b = l_s_t - l_s  # DUMMY

# Pedal travel
d_p = Q_("20 cm")  # DUMMY

# Stick deflection in handle (pitch)
d_e = Q_("0.1205 m")

# Stick deflection in handle (roll)
d_r = Q_("0.129 m")
