"""                  
Name: Wing
Department: Aerodynamics
Last updated: 05/06/2018 12:45 by Midas
"""

import os
import math as m
import numpy as np
import subprocess
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
FNULL = open(os.devnull,'w')

### Load distribution calculator
def computeloads():
    Datafile = "Airfoil"
    Re = 8000000
    M = 0.3128
    flapchordlength = 0.25
    maximumdeflectionangle = 30
    command_file = open("commands.in", "w")
    command_file.write('load ' + "../" + Datafile + ".dat\n\
    panel\n\
    panel\n\
    plop\n\
    g f\n\
    \n\
    gdes\n\
    f\n\
    " + str(1-flapchordlength) + "\n\
    999\n\
    0.5\n\
    " + str(maximumdeflectionangle) + "\n\
    x\n\
    \n\
    oper\n\
    visc " + str(Re) + "\n\
    M " + str(M) + "\n\
    type 1\n\
    pacc\n"\
    +  Datafile+"_results.dat\n\
    \n\
    iter\n 100 \n\
    seqp\n\
    aseq 0 25 0.5 \n\
    \n\
    \n\
    quit\n")
    command_file.close()

    run_xfoil_command = 'xfoil < ' + 'commands.in'
    subprocess.call(run_xfoil_command, stdout=FNULL, shell = True)

    data = np.genfromtxt(Datafile+'_results.dat',skip_header=12)
    alphaClmax = data[np.argmax(data[:,1]),0]

    subprocess.call('del ' + Datafile + '_results.dat', shell = True)

    print(alphaClmax)

    command_file = open("commands.in", "w")
    command_file.write('load ' + "../" + Datafile + ".dat\n\
    panel\n\
    panel\n\
    plop\n\
    g f\n\
    \n\
    gdes\n\
    f\n\
    " + str(1-flapchordlength) + "\n\
    999\n\
    0.5\n\
    " + str(maximumdeflectionangle) + "\n\
    x\n\
    \n\
    ppar\n\
    N\n\
    200\n\
    \n\
    \n\
    psav flappedairfoil.dat\n\
    oper\n\
    visc " + str(Re) + "\n\
    M " + str(M) + "\n\
    type 1\n\
    iter 1000 \n\
    seqp\n\
    a "  + str(alphaClmax) + "\n\
    cpwr cp.dat\n\
    \n\
    \n\
    quit\n")
    command_file.close()

    run_xfoil_command = 'xfoil < ' + 'commands.in'
    subprocess.call(run_xfoil_command, stdout= False, shell = True)

    pressures = np.genfromtxt('cp.dat',skip_header=3,skip_footer=1)

    Cy = 0
    Cx = 0
    Cm = 0
    xref = 0.25
    yref = 0.

    for i in range(np.size(pressures,0)-1):
        normal = np.array([pressures[i-1,2]-pressures[i+1,2],pressures[i+1,1]-pressures[i-1,1]])
        normalisednormal = normal/np.linalg.norm(normal)
        length = m.sqrt((pressures[i-1,1]-pressures[i+1,1])**2+(pressures[i-1,2]-pressures[i+1,2])**2)/2
        Cyi = normalisednormal[0]*pressures[i,0]*length
        Cxi = normalisednormal[1]*pressures[i,0]*length
        Cy += Cyi
        Cx += Cxi
        Cm += Cxi * (yref-pressures[i,2]) - Cyi * (xref-pressures[i,1])

    return Cy, Cx, Cm

print(computeloads())

#airfoil = np.genfromtxt(Datafile+'.dat')
#
# if not Keepresults[0].capitalize()=='Y':
#     os.system('del ' + Datafile+'_results.dat')
#
# alpha = data[:,0]
# Cl = data[:,1]
# Cd = data[:,2]
# Cdp = data[:,3]
# Cm = data[:,4]
# Top_Xtr = data[:,5]
# Bot_Xtr = data[:,6]
#
# command_file = open(file_path + Datafile + '_characteristics.dat', 'w')
# command_file.write('Maximum lift coefficient: ' + str(np.max(Cl)) + '\n\
# Minimum drag coefficient: ' + str(np.min(Cd)) + '\n\
# Maximum Cl/Cd: ' + str(np.max(Cl/Cd)) + '\n\
# Thickness: ' + str(np.max(airfoil[:,1])*2))
# command_file.close()