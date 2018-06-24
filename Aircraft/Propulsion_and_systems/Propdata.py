

import numpy as np
import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

data = np.loadtxt("../DataReal.txt", delimiter=" ", skiprows=1)


def data_reader(velocity, needed_value):
    if velocity.__class__ != int and velocity.__class__ != float:
        velocity.ito(ureg("m/s"))
        velocity = velocity.magnitude
    if velocity < 15:
        print("ERROR: please supply a value for velocity between 15 and 160 m/s")
        return
    index = int(velocity) - 15
    if needed_value == "Propeller efficiency":
        col = 0
        unit = "none"
    elif needed_value == "Total thrust":
        col = 1
        unit = "newton"
    elif needed_value == "Pitch angle":
        col = 2
        unit = "degree"
    elif needed_value == "Power needed":
        col = 3
        unit = "W"
    elif needed_value == "Torque":
        col = 4
        unit = "newton * meter"
    else:
        print("ERROR: please supply a valid column header as needed_value")
        return

    return Q_("{} {}".format(data[index, col], unit))
