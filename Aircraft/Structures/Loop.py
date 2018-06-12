'''This file is used to brute force those loops'''

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
import fileinput
import importlib
from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import numpy as np
from scipy import interpolate
import math as m
from Geometry import Geometry
from Geometry import Wing as GWing
import Wing
from Structures import Inertia
from Structures import Wing
from Structures import WingStress
from matplotlib import pyplot as plt



n = 10                      #number of the devided sections
b = Geometry.Wing.b/2         #Wing span
b = b.magnitude * ureg.meter
Normalstress = np.array([])
zarray = np.array([])

z = 0
while z < b.magnitude+0.1:
    NS = WingStress.Normal_stress_due_to_bending(0.15, Wing.airfoilordinate(0.15))
    Normalstress = np.append(Normalstress, NS.magnitude)
    zarray = np.append(zarray, z)
    text_to_search = 'z = ' + str(z)
    print(Wing.z, NS)

    z = z + b.magnitude/n

    replacement_text = 'z = ' + str(z)
    with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(text_to_search, replacement_text), end='')
    importlib.reload(Wing)
    importlib.reload(Inertia)
    importlib.reload(WingStress)

text_to_search = 'z = ' + str(z)
replacement_text = 'z = ' + str(0)
with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
    for line in file:
        print(line.replace(text_to_search, replacement_text), end='')

plt.plot(zarray, Normalstress)
plt.show()