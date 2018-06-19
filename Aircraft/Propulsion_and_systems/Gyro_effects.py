"""
Name: Gyroscopic effects
Department: Propulsion and Aircraft Systems
Last updated: 15/06/2018 13:58 by Ties
"""

"""
This file calculates gyroscopic effects of the propeller during turns and when applying a moment
"""

import sys
import scipy as sp
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

from Propulsion_and_systems import Propeller

# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centreline, pointing towards tail
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_rpm(inp):
    global rpm
    rpm = inp


# Start assigning values
rpm = Q_("2700 rpm")
rpm.ito(ureg("rad/min"))

# Make function which calculates gyro accelerations (prop only!) from moment and rate inputs


def input_moment(yaw_moment, pitch_moment, yaw_rate, pitch_rate):
    # Make compatible for inputs with and without units
    if yaw_moment.__class__ == int or yaw_moment.__class__ == float:
        yaw_moment *= ureg("newton*meter")
    if pitch_moment.__class__ == int or pitch_moment.__class__ == float:
        pitch_moment *= ureg("newton*meter")
    if yaw_rate.__class__ == int or yaw_rate.__class__ == float:
        yaw_rate *= ureg("rad/s")
    if pitch_rate.__class__ == int or pitch_rate.__class__ == float:
        pitch_rate *= ureg("rad/s")

    omega_y_dot = (yaw_moment - (Propeller.ixg - Propeller.izg) * pitch_rate * rpm) / Propeller.izg
    omega_z_dot = (pitch_moment - (Propeller.izg - Propeller.ixg) * yaw_rate * rpm) / Propeller.izg

    omega_y_dot.ito(ureg("rad/s**2"))
    omega_z_dot.ito(ureg("rad/s**2"))

    # print("Yaw acceleration: {}".format(omega_y_dot))
    # print("Pitch acceleration: {}".format(omega_z_dot))
    return [omega_y_dot, omega_z_dot]

# Test the function!
# input_moment(0, 0, Q_("5 deg/s"), 0)

# Make function which calculates gyro moments from rate and acceleration inputs


def input_acceleration(yaw_acceleration, pitch_acceleration, yaw_rate, pitch_rate):
    # Make compatible for inputs with and without units
    if yaw_acceleration .__class__ == int or yaw_acceleration.__class__ == float:
        yaw_acceleration *= ureg("rad/s**2")
    if pitch_acceleration .__class__ == int or pitch_acceleration.__class__ == float:
        pitch_acceleration *= ureg("rad/s**2")
    if yaw_rate.__class__ == int or yaw_rate.__class__ == float:
        yaw_rate *= ureg("rad/s")
    if pitch_rate.__class__ == int or pitch_rate.__class__ == float:
        pitch_rate *= ureg("rad/s")

    m_y = Propeller.izg * yaw_acceleration + (Propeller.ixg-Propeller.izg) * pitch_rate * rpm
    m_z = Propeller.izg * pitch_acceleration + (Propeller.izg - Propeller.ixg) * yaw_rate * rpm

    m_y.ito(ureg("newton * meter"))
    m_z.ito(ureg("newton * meter"))

    # print("Yaw moment: {}".format(m_y))
    # print("Pitch moment: {}".format(m_z))
    return [m_y, m_z]


input_acceleration(0, 0, 0, 1)
