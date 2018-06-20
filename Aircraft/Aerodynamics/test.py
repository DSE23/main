# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 17:32:23 2018

@author: Emma
"""
import numpy as np
import scipy as sp
import math
testset = np.array([(1,0),(10,1),(10**5,5),(10**15,15)])
x = np.log10((testset[:,0]))
y = (testset[:,1])
testtt = 100
curve = sp.interp(np.log10(testtt),y,x)
print('curve = ', 10**curve)


roughness_points = np.log10(np.array([(1.9e3,1e5),(6e5,1e4),(1e6,1.8e4),(7e6,9e4),(1e8,1.2e6),(9e8,1e7)])) #array for interpolation, from DATCOM graph
cutoffre = roughness_points[:,0] #x values selection
adm_rough = roughness_points[:,1] # y values selection
 #convert chord to inch for graph implementation
input_cut_re = np.log10(1e10)  #x-input for interpolation graph
curvere = 10**sp.interp(input_cut_re, adm_rough, cutoffre) #linear fit
print(curvere)

Re_used = 4e7
Cf_points = np.array([(3.5e5,0.0055),(5e5,0.005),(1.7e6,0.004),(9.5e6,0.003),(1.5e8,0.002),(1e9,0.0016),])
inputre = np.log10(Cf_points[:,0]) #xvals, logaritmhic scale
outputcf = Cf_points[:,1] #yvals
C_skinfric = (sp.interp(np.log10(Re_used), inputre,outputcf))
print(C_skinfric)
