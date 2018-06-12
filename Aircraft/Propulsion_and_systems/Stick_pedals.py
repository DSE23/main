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
    global l_s
    l_s = inp


def initialise_bottom_length(inp):
    global l_s_b
    l_s_b = inp


def initialise_pedal_travel(inp):
    global d_p
    d_p = inp


# End defining global variables

# Start assigning values to variables
# Stick length from hinge point to pilot's hand
l_s = Q_("81 cm")  # As defined by Gijs

# Stick length below hinge point
l_s_b = Q_("10 cm")  # DUMMY

# Pedal travel
d_p = Q_("30 cm")  # DUMMY
