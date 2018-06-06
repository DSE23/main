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
import scipy as sp
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

# Load case
# Max load factor (defined from requirements):
loadfactor = 20
# Gravitational acceleration
gravity = Q_("9.81 m/s**2")

# Determine max engine vertical force (in Newtons)
Engine.mass.ito(ureg.kg)
F_y_eng = Engine.mass * gravity * loadfactor
F_y_eng.ito(ureg.newton)

# Determine max prop vertical force (in Newtons)
Propeller.mass.ito(ureg.kg)
F_y_prop = Propeller.mass * gravity * loadfactor
F_y_prop.ito(ureg.newton)

print("F_y_eng = {}".format(F_y_eng))
print("F_y_prop = {}".format(F_y_prop))