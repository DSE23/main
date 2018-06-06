"""                  
Name: Wing
Department: Geometry
Last updated: 06/06/2018 12:38 by Boris 
"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

import numpy as np
import scipy as sp
import scipy.interpolate

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder

A = 5.5                         #Estimate aspect ratio
t = 0.4                         #Estimate taper
s = Q_("8.03 m")                #Estimate span (m)
Lambda25 = 0                    #Quarter chord sweep
CtoT = 0.15                     #Chord to thickness ratio
Spar2R = 0.18                     #Ratio of aileron at the root (postion of second spar)
Spar2T = 0.33                     #Ratio of aileron at the tip (postiion fo second spar)
Spar1R = 0.15                   #Ratio of LE HLD at root (postion of first spar)
Spar1T = 0.15                   #Ratio of LE HLD at tip (postion of first spar)
ChordR = Q_("2.015 m")         #Length of root (m)
ThSpar1 = Q_('5.0 mm')          #Thickness of Spar 1
ThSpar2 = Q_('5.0 mm')          #Thickness of Spar 2
ThSkin = Q_('3.0 mm')           #Thickness of the skin

airfoilcoordinates = np.genfromtxt("Airfoil.dat")    #Load coordinates
numberofcoordinates = np.size(airfoilcoordinates,0)  #Count total number of coordinates
airfoilinterpolant = sp.interpolate.interp1d(
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,0],
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,1],kind = 'cubic') #Interpolate

# Find the ordinate of the airfoil at an arbitrary position x, with 0 =< x =< 1
def airfoilordinate(x):
    return airfoilinterpolant(x)

TR = CtoR*ChordR                            #max thickness root in m
TT = TR*t                                   #max thickness tip in m


z = 0                                       #spanwise postion in meters
c = 0                                       #Chord wise postion in ratio

ChSpar1 = Spar1R + (Spar1T-Spar1R)*(z/s)    #Chord position of spar 1 with respect to
ChSpar2 = Spar2R + (Spar2T-Spar2R)*(z/s)    #Chord position of

## Here comes the function from Sam that relates chord to height

HSpar1= airfoilordinate(ChSpar1)*2
HSpar2= airfoilordinate(ChSpar2)*2

