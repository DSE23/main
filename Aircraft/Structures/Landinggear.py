'''
Simple landing gear sizing 
--------------------------
Boris Rowaan 
'''

import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from scipy import optimize
from scipy import interpolate
import math as m
from Geometry import Geometry

from Misc import ureg, Q_ # Imports the unit registry from the Misc folder

Chordlength = Q_('0.24 m')                                  #Lenght of chord landing gear
ThSkin = Q_('0.003 m')                                      #Thickness of the skin of the landing gear
lengthgear = Q_('0.8 m')
angle = 35                                                  #angle between vertical and landing gear
maxT = 0.217                                                #dimensionless max thickness of the airfoil on chordwise location


z = Q_('0 m')


#Material properties of the chosen material.
#Current chosen material:
#Epoxy/Carbon fiber, UD prepreg, QI lay-up
youngs_modulus = Q_("60.1 GPa")  #E
yield_strength = Q_("738 MPa")  #tensile
compr_strength = Q_("657 MPa") #compression
shear_modulus = Q_("23 GPa")   #G
poisson = 0.31                 # maximum 0.33
tau_max = Q_("35 MPa")
density = Q_("1560 kg/m**3")


'''Load geometry of the airfoil of the landing gear'''
##Ratio of height with respect to chord, airfoil coordinates
airfoilcoordinates = np.genfromtxt("../E475.dat",skip_header=1)    #Load coordinates
numberofcoordinates = np.size(airfoilcoordinates,0)  #Count total number of coordinates
airfoilinterpolant = sp.interpolate.interp1d(
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,0],
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,1],kind = 'cubic') #Interpolate

# Find the ordinate of the airfoil at an arbitrary position x, with 0 =< x =< 1
def airfoilordinate(x):
    return airfoilinterpolant(x)



def airfoilordinateroot(x,y):
    return airfoilordinate(x) - y

def airfoilinverse(y,location = 'TE'):
    airfoilrootinterpolant = sp.interpolate.interp1d(
        airfoilcoordinates[0:int(numberofcoordinates / 2) + 1, 0],
        airfoilcoordinates[0:int(numberofcoordinates / 2) + 1, 1]-y,kind = 'cubic')  # Interpolate
    if location == 'LE':
        x = 0.15-m.sqrt(0.15**2-y**2)
    elif location == 'TE':
        x = 0.999- 0.8/0.075*y
    return sp.optimize.newton(airfoilrootinterpolant,x)


'''----------------Centroid calculation---------------'''

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


def Area_Skin_x_c(Spar1, Spar2): #Input deminsionless chordwise location of spar 1 and spar 2
    n = 100 #number of sections
    dx = ((Spar2-Spar1)/n)
    x = Spar1
    Areaxc = 0
    arclength = 0
    for i in range(n):
        x = x + dx
        dxlength = dx * Chordlength
        dylength = abs(airfoilordinate(x - dx) - airfoilordinate(x)) * Chordlength
        dlength = np.sqrt(dxlength**2+dylength**2)
        dArea = dlength * ThSkin
        dAreaxc = dArea * x * Chordlength
        arclength += dlength
        Areaxc = Areaxc + dAreaxc
    Areaxc = Areaxc * 2                                         #Area times chord for both sides of the airfoil (therefore times 2)
    return Areaxc

AreaSkin = Area_Skin(0.001, 0.999)
AreaSkinxc = Area_Skin_x_c(0.001, 0.999)

centroid = (AreaSkinxc/AreaSkin)/Chordlength                    #dimensionless value of the centriod in chordlength (x)


'''-----------------Inertia calculations of the landing gear--------------------------------'''

def Calc_skin_inertia_Ixx(Spar1, Spar2):
    n = 100  # number of sections
    dx = ((Spar2 - Spar1) / n)
    x = Spar1
    Ixx = 0
    for i in range(n):
        x = x + dx
        dxlength = dx * Chordlength
        y = ((airfoilordinate(x - dx) + airfoilordinate(x)) / 2) * Chordlength
        dy = abs(airfoilordinate(x - dx) - airfoilordinate(x)) * Chordlength
        length = np.sqrt(dxlength ** 2 + dy ** 2)
        dIxx = length * ThSkin * (y ** 2)
        Ixx = Ixx + dIxx
    Ixx = Ixx * 2
    return Ixx

def Calc_skin_inertia_Iyy(Spar1, Spar2):
    n = 100  # number of sections
    dx = ((Spar2 - Spar1) / n)
    x = Spar1
    Iyy = 0
    for i in range(n):
        x = x + dx
        xlength = x * Chordlength
        dxlength = dx * Chordlength
        y = ((airfoilordinate(x - dx) + airfoilordinate(x)) / 2) * Chordlength
        dy = abs(airfoilordinate(x - dx) - airfoilordinate(x)) * Chordlength
        length = np.sqrt(dxlength ** 2 + dy ** 2)
        dIyy = length * ThSkin * ((abs(xlength - (centroid * Chordlength))) ** 2)
        Iyy = Iyy + dIyy
    Iyy = Iyy * 2
    return Iyy

IxxSkin = Calc_skin_inertia_Ixx(0.001, 0.999)
IyySkin = Calc_skin_inertia_Iyy(0.001, 0.999)
IxySkin = Q_('0 m**4')


def landing_stress(z):

    '''------------------Loading calculations------------------------'''
    load = Geometry.Masses.W_MTOW * Q_('9.81 m / s**2') / 2                         #Force exerted per landing gear
    G = 10                                                           #Load factor
    G_load = load * G
    B_force = np.sin(np.radians(angle)) * G_load                           #bending force
    C_force = np.cos(np.radians(angle)) * G_load                           #compression force

    Mx = (lengthgear - z) * B_force
    My = Q_('0 N * m')

    '''-------------------Stress calculations-------------------------'''

    def Normal_stress_due_to_bending(cs, y): # Normal stress due to bending
        denominator_inertia_term = IxxSkin*IyySkin-IxySkin**2
        inertia_term_1 = (IyySkin*y-IxySkin*cs)/denominator_inertia_term
        inertia_term_2 = (IxxSkin*cs-IxySkin*y)/denominator_inertia_term
        sigma_zs = My*inertia_term_1 + Mx*inertia_term_2
        strain = sigma_zs / youngs_modulus
        return sigma_zs, strain #Gives the normal stress function for a given span zs, and x- and y- coordinate

    sigma, strain = Normal_stress_due_to_bending(maxT, airfoilordinate(maxT))
    strain.ito(ureg('m**-1'))
    return sigma, strain
'''-----------------Stress calculation loop-----------------------------'''

sigmalist = np.array([])                            #defining some empty lists
strainlist = np.array([])
zlist = np.array([])
n = 40                                              #number of sections

while z <= lengthgear:
    sigma, strain = landing_stress(z)
    zlist = np.append(zlist, z)
    sigmalist = np.append(sigmalist, sigma)
    strainlist = np.append(strainlist, strain)
    z = z + (lengthgear/n)

Vol = AreaSkin * Chordlength
Mass = Vol * density
print(Mass)

plt.plot(zlist, strainlist)
plt.show()
