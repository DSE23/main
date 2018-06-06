"""
Name: Engine Mounts
Department: Propulsion and Aircraft Systems
Last updated: 05/06/2018 16:52 by Ties
"""

"""
This file calculates the maximum forces experienced by the engine mounts, 
for use by the structures department
"""

import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

# Import data
import Geometry
import Engine
import Propeller

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, backwards
# Y axis: up
# Z axis: left

