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
from Structures import Inertia
from Structures import Wing
from Structures import WingStress
from Structures import Shear
from matplotlib import pyplot as plt



n = 10                      #number of the devided sections
b = Wing.s         #Wing span
b = b.magnitude * ureg.meter
Normalstress = np.array([])

Vol_mat_spar1 = Q_('0 m**3')
Vol_mat_spar2 = Q_('0 m**3')
Vol_mat_skin = Q_('0 m**3')
Vol_mat_string = Q_('0 m**3')
Vol_mat_clamp = Q_('0 m**3')
Vol_mat_wing = Q_('0 m**3')

Old_N_stringers = Wing.N_stringers

Lmomentlist = np.array([])
Ixxlist = np.array([])
Iyylist = np.array([])
zarray = np.array([])
Farray = np.array([])

z = 0
while z <= b.magnitude:
    NS = WingStress.Normal_stress_due_to_bending(0.18, Wing.airfoilordinate(0.18))[0]
    Normalstress = np.append(Normalstress, NS.magnitude)
    zarray = np.append(zarray, z)
    F = Shear.F
    Farray = np.append(Farray, F)
    z *= Q_('m')
    L_moment = WingStress.computeloads(z)[3]
    print('L_moment', L_moment)
    z = z.magnitude
    '''Calculate all the subweights of the wing '''
    Vol_mat_spar1 = Vol_mat_spar1 + Wing.AreaSpar1*(b/n)
    Vol_mat_spar2 = Vol_mat_spar2 + Wing.AreaSpar2 * (b / n)
    Vol_mat_skin = Vol_mat_skin + Wing.Area_Skin(Wing.ChSpar1, Wing.ChSpar2) * (b / n)
    Vol_mat_string = Vol_mat_string + Wing.AreaStringers * (b / n)
    Vol_mat_wing = Vol_mat_wing + Wing.Area * (b / n)
    Lmomentlist = np.append(Lmomentlist, L_moment)
    Ixxlist = np.append(Ixxlist, Inertia.Ixx_wb)
    Iyylist = np.append(Iyylist, Inertia.Iyy_wb)
    Dist_between_spars = Wing.ChSpar2*Wing.Chordlength - Wing.ChSpar1*Wing.Chordlength

    print(L_moment, Dist_between_spars)
    print(Wing.z, NS, Wing.N_stringers)

    text_to_search = 'z = ' + str(z)
    z = z + b.magnitude/n
    replacement_text = 'z = ' + str(z)
    with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(text_to_search, replacement_text), end='')

    if Geometry.Fuselage.b_f.magnitude < z < 1.4:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 2
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    if 1.4 < z < 3.0:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 20
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    if 3.0 < z < 3.5:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 2
        replacement_text = 'N_stringers = ' + str(New_N_stringers)
        with fileinput.FileInput('Wing.py', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

    if z > 3.5:
        text_to_search = 'N_stringers = ' + str(Wing.N_stringers)
        New_N_stringers = 2
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

Density = Q_('1560 kg / m**3')


Weightspar1 = Density *Vol_mat_spar1
Weightspar2 = Density * Vol_mat_spar2
Weightskin = Density * Vol_mat_skin
Weightstring = Density * Vol_mat_string
Weightwing = Density * Vol_mat_wing

print("Wingweight half span main wing", Weightwing)

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

print('debug:', len(zarray), len(Iyylist))

# plot with various axes scales
fig = plt.figure()

# linear
plt.subplot(2, 2, 1)
plt.plot(zarray, Ixxlist)

plt.subplot(2, 2, 2)
plt.plot(zarray, Farray)

plt.subplot(2, 2, 3)
plt.plot(zarray, Lmomentlist)

plt.subplot(2, 2, 4)
plt.plot(zarray, Normalstress)

plt.show()

