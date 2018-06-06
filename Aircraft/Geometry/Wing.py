"""                  
Name: Wing
Department: Geometry
Last updated: 06/06/2018 12:38 by Boris 
"""

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