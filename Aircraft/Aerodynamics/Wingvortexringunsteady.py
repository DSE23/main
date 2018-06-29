"""
Name: Airfoil discrete vortex simulation
Department: Aerodynamics
Last updated: 11/06/2018 12:45 by Sam
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import time

N = 40

### Input geometry
fuselagediameter = 1.04
wingspan = 7.54
rootchord = 1.8
tipchord = 0.8



### Output geometry
fuselagehalfidameter = fuselagediameter/2
winghalfspan = wingspan/2
effectivewinghalfspan = winghalfspan-fuselagehalfidameter
winggeometriccenter = fuselagehalfidameter+effectivewinghalfspan/2

theta = np.linspace(0,m.pi,N,endpoint=True)
x = np.linspace(0,1,40)
y = winggeometriccenter - effectivewinghalfspan*np.cos(theta)/2

Y,X = np.meshgrid(y,x)

def airfoilcamberline(x):
    return 0*x

def chordlength(y):
    return rootchord - (rootchord-tipchord)/winghalfspan*y

def chordscaling(x,y):
    return (x-0.25)*chordlength(y)

panelarray = np.zeros((39*39,8,3))
panelarray[:,0,0] = X[0:39,0:39].reshape(-1)
panelarray[:,0,1] = Y[0:39,0:39].reshape(-1)
panelarray[:,1,0] = X[1:40,1:40].reshape(-1)
panelarray[:,1,1] = Y[0:39,0:39].reshape(-1)
panelarray[:,2,0] = X[1:40,1:40].reshape(-1)
panelarray[:,2,1] = Y[1:40,1:40].reshape(-1)
panelarray[:,3,0] = X[0:39,0:39].reshape(-1)
panelarray[:,3,1] = Y[1:40,1:40].reshape(-1)
panelarray[:,0:4,0] = chordscaling(panelarray[:,0:4,0],panelarray[:,0:4,1])
panelarray[:,0:4,2] = airfoilcamberline(panelarray[:,0:4,0])
panelarray[:,4,:] = panelarray[:,2,:] - panelarray[:,0,:]
panelarray[:,5,:] = panelarray[:,3,:] - panelarray[:,1,:]
panelarray[:,6,:] = np.cross(panelarray[:,4,:],panelarray[:,5,:])
panelarray[:,7,:] = panelarray[:,6,:]//(np.sqrt((panelarray[:,6,:]*panelarray[:,6,:]).sum(axis=1)).reshape(-1,1))

