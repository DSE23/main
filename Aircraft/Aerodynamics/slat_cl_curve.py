# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 15:32:08 2018

@author: Emma
"""

import numpy as np
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
import matplotlib.pyplot as plt

slope = np.genfromtxt('../MS18_for_slat.txt')
x = slope[:,0]
y = slope[:,1]

x_slat = slope[18:40,0] + 4

y_slat = slope[18:40,1] + 0.531


plt.plot(x,y,'k')
plt.plot(x_slat,y_slat,'k--')
plt.title('Lift curve slope with slats deployed')
plt.xlabel(r'$\alpha$ [°]')
plt.ylabel(r'$C_{l}$ [-]')
plt.legend(['Clean', 'Slats deployed'])
plt.show