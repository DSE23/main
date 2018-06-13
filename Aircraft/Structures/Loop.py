'''This file is used to brute force those loops'''

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
import fileinput
import importlib
from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder
import numpy as np
from scipy import interpolate
import math as m
import Wing
from Geometry import Geometry
from Geometry import Wing as GWing
from Structures import Inertia
from Structures import Wing
from Structures import WingStress
from matplotlib import pyplot as plt



n = 10                      #number of the devided sections
b = Geometry.Wing.b/2         #Wing span
b = b.magnitude * ureg.meter
Normalstress = np.array([])
zarray = np.array([])
Vol_mat_spar1 = Q_('0 m**3')
Vol_mat_spar2 = Q_('0 m**3')
Vol_mat_skin = Q_('0 m**3')
Vol_mat_string = Q_('0 m**3')
Vol_mat_wing = Q_('0 m**3')

Old_N_stringers = Wing.N_stringers

Lmomentlist = np.array([])

z = 0
while z < b.magnitude+0.1:
    NS = WingStress.Normal_stress_due_to_bending(0.15, Wing.airfoilordinate(0.15))
    Normalstress = np.append(Normalstress, NS.magnitude)
    zarray = np.append(zarray, z)

    '''Calculate all the subweights of the wing '''
    Vol_mat_spar1 = Vol_mat_spar1 + Wing.AreaSpar1*(b/n)
    Vol_mat_spar2 = Vol_mat_spar2 + Wing.AreaSpar2 * (b / n)
    Vol_mat_skin = Vol_mat_skin + Wing.Area_Skin(Wing.ChSpar1, Wing.ChSpar2) * (b / n)
    Vol_mat_string = Vol_mat_string + Wing.AreaStringers * (b / n)
    Vol_mat_wing = Vol_mat_wing + Wing.Area * (b / n)
    Lmomentlist = np.append(Lmomentlist, WingStress.L_moment)
    print(WingStress.L_moment)
    print(Wing.z, NS, Wing.N_stringers)
    text_to_search = 'z = ' + str(z)
    z = z + b.magnitude/n
    replacement_text = 'z = ' + str(z)
    with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(text_to_search, replacement_text), end='')

    if Geometry.Fuselage.b_f.magnitude < z < 2.0:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 10
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    if 2.0 < z < 2.5:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 10
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    if 2.5 < z < 3.5:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 10
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    if z > 3.5:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 10
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    importlib.reload(Wing)
    importlib.reload(Inertia)
    importlib.reload(WingStress)

text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
replacement_text = 'N_stringers = ' + str(Old_N_stringers)
with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
    for line in file:
        print(line.replace(text_to_search, replacement_text), end='')

text_to_search = 'z = ' + str(z)
replacement_text = 'z = ' + str(0)
with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
    for line in file:
        print(line.replace(text_to_search, replacement_text), end='')

'''Random Density variable here'''

Density = Q_('1700 kg / m**3')


Weightspar1 = Density *Vol_mat_spar1
Weightspar2 = Density * Vol_mat_spar2
Weightskin = Density * Vol_mat_skin
Weightstring = Density * Vol_mat_string
Weightwing = Density * Vol_mat_wing



# with is like your try .. finally block in this case
with open('StrucVal.py', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

data[6] = 'Weightspar1 = Q_(\"' + str(Weightspar1) + '\")\n'
data[7] = 'Weightspar2 = Q_(\"' + str(Weightspar2) + '\")\n'
data[8] = 'Weightskin = Q_(\"' + str(Weightskin) + '\")\n'
data[9] = 'Weightstring = Q_(\"' + str(Weightstring) + '\")\n'
data[10] = 'Weightwing = Q_(\"' + str(Weightwing) + '\")\n'
data[17] = 'Density = Q_(\"' + str(Density) + '\")\n'

# and write everything back
with open('StrucVal.py', 'w') as file:
    file.writelines(data)

plt.plot(zarray, Normalstress)
plt.show()