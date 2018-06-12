"""
Name: Linkages
Department: Propulsion and Aircraft Systems
Last updated: 12/06/2018 11:38 by Ties
"""

"""
This file calculates the forces and moments on the control surfaces, linkages and stick 
"""

import sys
import scipy as sp
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

from Geometry import Geometry
import numpy as np

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_max_roll_stick_force(inp):
    global f_s_r_max
    f_s_r_max = inp


def initialise_stick_length(inp):
    global l_stick
    l_stick = inp


def initialise_rudder_travel(inp):
    global l_rudder
    l_rudder = inp


# End defining global variables

# Start assigning values to variables
# Maximum stick force in roll
f_s_r_max = Q_("20 lbf")  # From MIL-F-8785c

# Maximum stick force in pitch
f_s_p_max = f_s_r_max*2  # From 1:2:3 ratio for control forces

# Maximum rudder force
f_r_max = f_s_r_max*3  # From 1:2:3 ratio for control forces

# Stick length
l_stick = Q_("81 cm")  # As defined by Gijs

# Rudder travel
l_rudder = Q_("30 cm")  # DUMMY
# End assigning values


# Calculate maximum hinge moment in roll
m_s_r_max = f_s_r_max * l_stick
m_s_r_max.ito(ureg("newton * meter"))
print(m_s_r_max)

# Calculate maximum hinge moment in pitch
m_s_p_max = f_s_p_max * l_stick
m_s_p_max.ito(ureg("newton * meter"))
print(m_s_p_max)

# Calculate maximum hinge moment in yaw
r_rudder = l_rudder/np.tan(Geometry.V_tail.delta_r)
m_r_max = f_r_max * r_rudder
m_r_max.ito(ureg("newton * meter"))
print(m_r_max)
