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

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_stick_moment(inp):
    global m_s
    m_s = inp


# End defining global variables

# Start assigning values to variables
m_s = Q_("86 N")
