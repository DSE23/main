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
from math import *

from Misc import ureg, Q_ # Imports the unit registry fron the Misc folder

A = 5.5                         #Estimate aspect ratio
t = 0.4                         #Estimate taper
s = Q_("8.03 m")                #Estimate span (m)
Lambda25 = 0                    #Quarter chord sweep
CtoT = 0.15                     #Max Chord to thickness ratio
Spar2R = 1-0.18                      #Chordwise location of second spar at the root
Spar2T = 1-0.33                     #Chordwise location of second spar at the tip
Spar1R = 0.15                   #Chordwise location of first spar at the root
Spar1T = 0.15                   #Chordwise location of first spar at the tip
ChordR = Q_("2.015 m")         #Length of root (m)
ThSpar1 = Q_('0.005 m')          #Thickness of Spar 1
ThSpar2 = Q_('0.005 m')          #Thickness of Spar 2
ThSkin = Q_('0.003 m')           #Thickness of the skin

##Assumption
centroid = 0.5

##Stringers                     # C stringer dimentions
h_str = Q_('0.025 m')            # height of the stringer
w_str = Q_('0.025 m')            #width of the stringer
t_str = Q_('3.0 mm')            #thickness of the stringer


z = Q_('0 m')                                       #spanwise posotion in meters
c = 0                                               #Chord wise postion in ratio



##Ratio of height with respect to chord, airfoil coordinates
airfoilcoordinates = np.genfromtxt("../Airfoil.dat")    #Load coordinates
numberofcoordinates = np.size(airfoilcoordinates,0)  #Count total number of coordinates
airfoilinterpolant = sp.interpolate.interp1d(
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,0],
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,1],kind = 'cubic') #Interpolate

# Find the ordinate of the airfoil at an arbitrary position x, with 0 =< x =< 1
def airfoilordinate(x):
    return airfoilinterpolant(x)

TR = airfoilordinate(0.15)*ChordR                            #max thickness root in m
TT = TR*t                                   #max thickness tip in m

def Chord_loc_Spar(zs,SparR,SparT):             #input spanwise location in m and
    ChSpar = SparR + (SparT-SparR)*(zs/s)  #Chord position of spar 1 with respect to leading edge
    return ChSpar


ChSpar1 = Chord_loc_Spar(z,Spar1R,Spar1T)
ChSpar2 = Chord_loc_Spar(z,Spar2R,Spar2T)

# print(ChSpar1)
# print(ChSpar2)

def length_chord(zs):             #length of chord with respect to the spanwise postion
    lengthchord = ChordR*(1-((1-t)*(zs/s)))
    return lengthchord
Chordlength = length_chord(z)


## Here comes the function from Sam that relates chord to height

def H_in_m(x,zs):                             #input chord and span respectively
    HSpar= airfoilordinate(x)*ChordR*(1-(1-t)*(zs/s))         #Function to calculate the height in m measured from the line symmetry depending on the chord and span
    return HSpar

HSpar1 = H_in_m(ChSpar1, z)*2        #height of spar 1
HSpar2 = H_in_m(ChSpar2, z)*2        #height of spar 2

# print(HSpar1)
# print(HSpar2)

def Angle(cs):                          #input chord ratio
    n = 100  # number of sections
    dx = 1/n
    dy=airfoilordinate(cs+dx)-airfoilordinate(cs)
    angle = tan(dy/dx)
    angle *= Q_('rad')
    return angle

## Calculation of centriod

#Area of Spars
AreaSpar1 = HSpar1 * ThSpar1
AreaSpar2 = HSpar2 * ThSpar2

#Area of Skin
def Area_Skin(Spar1, Spar2):                            #Input deminsionless chordwise location of spar 1 and spar 2
    n = 100 #number of sections
    dx = ((Spar2-Spar1)/n)
    x = Spar1
    Area = 0
    for i in range(n):
        x = x + dx
        dxlength = dx * Chordlength
        dylength = abs(airfoilordinate(x - dx) - airfoilordinate(x)) * Chordlength
        dlength = np.sqrt(dxlength**2+dylength**2)
        dArea = dlength * ThSkin
        Area = Area + dArea
    Area = Area * 2                                 # Area of both sides of the airfoil
    return Area

#Area of the Stringers


## Area multiplied with the distance from the reference point (leading edge c=0)

#For the spars
AreaSpar1xc = AreaSpar1 * ChSpar1 * Chordlength
AreaSpar2xc = AreaSpar2 * ChSpar2 * Chordlength

def Area_Skin_x_c(Spar1, Spar2):                            #Input deminsionless chordwise location of spar 1 and spar 2
    n = 100 #number of sections
    dx = ((Spar2-Spar1)/n)
    x = Spar1
    Areaxc = 0
    for i in range(n):
        x = x + dx
        dxlength = dx * Chordlength
        dylength = abs(airfoilordinate(x - dx) - airfoilordinate(x)) * Chordlength
        dlength = np.sqrt(dxlength**2+dylength**2)
        dArea = dlength * ThSkin
        dAreaxc = dArea * x * Chordlength
        Areaxc = Areaxc + dAreaxc
    Areaxc = Areaxc * 2                                         #Area times chord for both sides of the airfoil (therefore times 2)
    return Areaxc


def area_stringers(n_st):
    A_1 = h_str*t_str
    A_2 = (w_str-t_str)*t_str
    A_stringer = A_1 + A_2
    A_tot_stringer = A_stinger * n_st
    return A_tot_stringer

def stif_loc(z, t_sk, n_st, x):
    total_perimeter = sp.integrate.quad(Wing.airfoilordinate(x), Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), Wing.Chord_loc_Spar(z,Wing.Spar2R,Wing.Spar2T)) #m

    spacing = total_perimeter / ((n_st + 1) / 2)
    x_y_angle_coords = []
    for i in range(n_st):
        local_spacing = i * spacing
        rot_angle = Wing.Angle(x) + radians(180)
        x_coordinate = Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T) + sp.integrate.quad(math.cos(Wing.angle(x)),Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), localspacing) 
        #x_coordinate = (-1) * (local_spacing - circle_perim) * cos(atan(0.5 * h / (C_a - 0.5 * h)))
        y_coordinate = Wing.airfoilordinate(Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T)+local_spacing)-sp.integrate.quad(math.sin(Wing.angle(x)),Wing.Chord_loc_Spar(z,Wing.Spar1R,Wing.Spar1T), localspacing)
        
        apnd_itm = (x_coordinate, y_coordinate, rot_angle)
        x_y_angle_coords.append(apnd_itm)
        apnd_itm = (x_coordinate, -y_coordinate, -rot_angle)
        x_y_angle_coords.append(apnd_itm)

        # print "Stif.", i, "\t x:", x_coordinate, "\t y:", y_coordinate, "\t angle:", degrees(rot_angle)

    return x_y_angle_coords  # [(stringer0 x,y,rot),(stringer1 x,y,rot), ...]


Print("je moeder") 





