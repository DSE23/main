import numpy as np

array = np.genfromtxt("dataforjordy.dat")

cn = array[:,1]*np.cos(np.radians(array[:,0]))+array[:,2]*np.cos(np.radians(array[:,0]))
cp = 0.25- array[:,4]/cn

print(np.max(cp),np.min(cp))

print(cp)