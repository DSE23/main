"""                  
Name: Wing
Department: Structures
Last updated: 05/06/2018 12:45 by Boris
"""
# Midas is een held
# Import packages

from pint import UnitRegistry
ureg = UnitRegistry()

## Defining variables

A = 5.5                         #Estimate aspect ratio
t = 0.4                         #Estimate taper
s = 8.03                        #Estimate span (m)
s *= ureg.meter                 #giving span meter as unit
Lambda25 = 0                    #Quarter chord sweep
CtoT = 0.15                     #Chord to thickness ratio
Spar2R = 0.18                     #Ratio of aileron at the root (postion of second spar)
Spar2T = 0.33                     #Ratio of aileron at the tip (postiion fo second spar)
Spar1R= 0.15                   #Ratio of LE HLD at root (postion of first spar)
Spar1T= 0.15                   #Ratio of LE HLD at tip (postion of first spar)
ChordR = 2.015                  #Length of root (m)
ChordR *= ureg.meter            #giving root chord meter as unit

z = 0                           #spanwise postion in meters

TR = CtoR*ChordR                #max thickness root in m
TT = TR*t                       #max thickness tip in m
ChSpar1 = Spar1R + (Spar1T-Spar1R)*(z/s)   #Chord position of spar 1 with respect to
ChSpar2 = Spar2R + (Spar2T-Spar2R)*(z/s)   #Chord position of
CtoT =
HSpar1= CtoT - (ChSpar1-CposmaxT)
HSpar2=






print(Lambda25)
