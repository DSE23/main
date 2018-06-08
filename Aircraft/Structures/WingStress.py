"""                  
Name: Wing_Stress_Calculations 
Department: Structures 
Last updated: 08/06/2018 12:38 by Boris 
"""

## Forces
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

import numpy as np
import scipy as sp
from scipy import interpolate
import math as m
import Geometry
import Structures
from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder


cl = 1
cd = 0.05
cm = 0.05
b = Geometry.Wing.b         #Wing span
z = Structures.Wing.z      #Span wise postion of the wing in m
ChordR = Geometry.Wing.c_r      #root chord in m

av_chord = (Wing.length_chord(z)+Wing.length_chord(b))/2        #average chord right from the crossection (m)
spanleft = b - z
print(av_chord)
print(spanleft)