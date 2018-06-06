"""                  
Name: Wing
Department: Structures
Last updated: 05/06/2018 12:45 by Boris
"""
# Midas is een held
# Import packages

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder

## Defining variables



z = 0                                       #spanwise postion in meters

TR = CtoR*ChordR                            #max thickness root in m
TT = TR*t                                   #max thickness tip in m
ChSpar1 = Spar1R + (Spar1T-Spar1R)*(z/s)    #Chord position of spar 1 with respect to
ChSpar2 = Spar2R + (Spar2T-Spar2R)*(z/s)    #Chord position of
CtoT =
HSpar1= CtoT - (ChSpar1-CposmaxT)
HSpar2=

print(Lambda25)
