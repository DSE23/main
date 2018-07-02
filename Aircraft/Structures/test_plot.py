# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 18:27:59 2018

@author: tobia
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D

cs = np.genfromtxt("cs_data.dat")
zarray = np.genfromtxt("zarray_data.dat")
clist = np.genfromtxt("y_data.dat")
hlist = np.genfromtxt("z_data.dat")
Farray = np.genfromtxt("Farray_data.dat")


x = zarray
y = clist
z = hlist


cm = plt.get_cmap('rainbow')
if min(cs) >= max(cs):
    cNorm = matplotlib.colors.Normalize(vmin=-min(cs), vmax=min(cs))
if min(cs) < max(cs):
    cNorm = matplotlib.colors.Normalize(vmin=-max(cs), vmax=max(cs))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
scalarMap.set_array(cs)
fig.colorbar(scalarMap)
plt.show()

F = Farray

cm = plt.get_cmap('plasma_r')
if min(F) >= max(F):
    cNorm = matplotlib.colors.Normalize(vmin=min(F), vmax=max(F))
if min(F) < max(F):
    cNorm = matplotlib.colors.Normalize(vmin=min(F), vmax=max(F))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z, c=scalarMap.to_rgba(F))
scalarMap.set_array(F)
fig.colorbar(scalarMap)
plt.show()