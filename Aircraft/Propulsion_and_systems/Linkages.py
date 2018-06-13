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
import Stick_pedals
from Structures import StrucVal

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


def initialise_elevator_pushrod_diameter(inp):
    global d_pushrod_e
    d_pushrod_e = inp


def initialise_elevator_pushrod_length(inp):
    global l_pushrod_e
    l_pushrod_e = inp


# End defining global variables

# Start assigning values to variables
# Maximum stick force in roll
f_r_max = Q_("20 lbf")  # From MIL-F-8785c

# Maximum stick force in pitch
f_p_max = f_r_max*2  # From 1:2:3 ratio for control forces

# Maximum yaw force on pedals
f_y_max = f_r_max*3  # From 1:2:3 ratio for control forces

# End assigning values

# Hinge moments calculation

# Calculate maximum hinge moment in roll
m_r_max = f_r_max * Stick_pedals.l_s
m_r_max.ito(ureg("newton * meter"))
print("Maximum roll hinge moment: {}".format(m_r_max))

# Calculate maximum hinge moment in pitch
m_p_max = f_p_max * Stick_pedals.l_s
m_p_max.ito(ureg("newton * meter"))
print("Maximum pitch hinge moment: {}".format(m_p_max))

# Arm between rudder and cable/pushrod attachment
r_rudder = Stick_pedals.d_p/np.tan(Geometry.V_tail.delta_r)
# Calculate maximum hinge moment in yaw
m_y_max = f_y_max * r_rudder
m_y_max.ito(ureg("newton * meter"))
print("Maximum yaw hinge moment: {}".format(m_y_max))

# Calculate forces in pushrods
# Max roll force on stick in pushrod attachment
f_r_b_max = (Stick_pedals.l_s/Stick_pedals.l_s_b) * f_r_max
f_r_b_max.ito(ureg("newton"))
print("Maximum roll force on stick in pushrod attachment: {}".format(f_r_b_max))

# Max pitch force on stick in pushrod attachment
f_p_b_max = (Stick_pedals.l_s/Stick_pedals.l_s_b) * f_p_max
f_p_b_max.ito(ureg("newton"))
print("Maximum pitch force on stick in pushrod attachment: {}".format(f_p_b_max))

# Max yaw force in cables equals f_y_max above
f_y_max.ito(ureg("newton"))
print("Maximum yaw force on rudder cables: {}".format(f_y_max))

# Buckling calculation of pushrods
# Define dimensions
# Elevator pushrod diameter
d_pushrod_e = Q_("10 cm")
# Elevator pushrod length
l_pushrod_e = Q_("3 m")

# Calculate critical moment of inertia for buckling
i_xx_e = (f_p_max*l_pushrod_e)/(np.pi**2*StrucVal.youngs_modulus)

# Calculate skin thickness needed for critical moment of inertia
t_pushrod_e = d_pushrod_e/2 - (i_xx_e*4/np.pi - (d_pushrod_e/2)**4)**(1/4)

print("Minimum skin thickness needed to prevent buckling at {} diameter: {}".format(d_pushrod_e,t_pushrod_e))
