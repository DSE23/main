"""                  
Name: Wing
Department: Geometry
Last updated: 06/06/2018 12:38 by Boris 
"""
# x is chordwise from LE to TE. y is positive upwards.
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

## MAIN WING

A = Geometry.Wing.A                         #Estimate aspect ratio
t = Geometry.Wing.taper                         #Estimate taper
s = Geometry.Wing.b/2-Geometry.Wing.horn                #Estimate span (m)
Lambda25 = 0                    #Quarter chord sweep
CtoT = 0.15                     #Max Chord to thickness ratio
Spar2R = (1-Geometry.Wing.c_a/Geometry.Wing.c_r)  #Chordwise location of second spar at the root
Spar2T = (1-Geometry.Wing.c_a/Geometry.Wing.c_t)  #Chordwise location of second spar at the tip
Spar1R = 0.18                   #Chordwise location of first spar at the root
Spar1T = 0.18                   #Chordwise location of first spar at the tip
ChordR = Geometry.Wing.c_r      #Length of root (m)
ThSpar1 = Q_('0.002 m')         #Thickness of Spar 1
ThSpar2 = Q_('0.002 m')        #Thickness of Spar 2
ThSkin = Q_('0.0045 m')         #Thickness of the skin
N_stringers = 2                 #Number of stringers
ClampH = Q_('0.01 m')           #height of the clamps at the top of the spars
ClampW = Q_('0.01 m')           #width of the clamps at the top of the spars


##Stringers                     # C stringer dimentions
h_str = Q_('0.025 m')            # height of the stringer
w_str = Q_('0.025 m')            #width of the stringer
t_str = Q_('0.003 m')            #thickness of the stringer


z = 0                       #spanwise posotion in meters
z *= Q_('meter')
c = 0                                               #Chord wise postion in ratio

# ## VERTICAL TAIL
# 
#VTA = Geometry.Wing.A                         #Estimate aspect ratio
#VTt = Geometry.Wing.taper                         #Estimate taper
#VTs = Geometry.Wing.b/2-Geometry.Wing.horn                #Estimate span (m)
#VTLambda25 = 0                    #Quarter chord sweep
#VTCtoT = 0.15                     #Max Chord to thickness ratio
#VTSpar2R = (1-Geometry.H_tail.cr_r)  #Chordwise location of second spar at the root
#VTSpar2T = (1-Geometry.H_tail.cr_r)  #Chordwise location of second spar at the tip
#VTSpar1R = 0.18                   #Chordwise location of first spar at the root
#VTSpar1T = 0.18                   #Chordwise location of first spar at the tip
#VTChordR = Geometry.Wing.c_r      #Length of root (m)
#VTThSpar1 = Q_('0.002 m')         #Thickness of Spar 1
#VTThSpar2 = Q_('0.002 m')        #Thickness of Spar 2
#VTThSkin = Q_('0.0018 m')         #Thickness of the skin
#VTN_stringers = 2                 #Number of stringers
#VTClampH = Q_('0.03 m')           #height of the clamps at the top of the spars
#VTClampW = Q_('0.03 m')           #width of the clamps at the top of the spars
#
#
###Stringers                     # C stringer dimentions
#VTh_str = Q_('0.025 m')            # height of the stringer
#VTw_str = Q_('0.025 m')            #width of the stringer
#VTt_str = Q_('0.003 m')            #thickness of the stringer
#
# 
#
# ## HORIZONTAL TAIL
# 
#HTA = Geometry.H_tail.A                         #Estimate aspect ratio
#HTt = Geometry.H_tail.taper                         #Estimate taper
#HTs = Geometry.H_tail.b/2-Geometry.Wing.horn                #Estimate span (m)
#HTLambda25 = 0                    #Quarter chord sweep
#HTCtoT = 0.15                     #Max Chord to thickness ratio
#HTSpar2R = (1-Geometry.H_tail.ce_c)  #Chordwise location of second spar at the root
#HTSpar2T = (1-Geometry.H_tail.ce_c)  #Chordwise location of second spar at the tip
#HTSpar1R = 0.18                   #Chordwise location of first spar at the root
#HTSpar1T = 0.18                   #Chordwise location of first spar at the tip
#HTChordR = Geometry.H_tail.c_r      #Length of root (m)
#HTThSpar1 = Q_('0.002 m')         #Thickness of Spar 1
#HTThSpar2 = Q_('0.002 m')        #Thickness of Spar 2
#HTThSkin = Q_('0.0018 m')         #Thickness of the skin
#HTN_stringers = 2                 #Number of stringers
#HTClampH = Q_('0.03 m')           #height of the clamps at the top of the spars
#HTClampW = Q_('0.03 m')           #width of the clamps at the top of the spars
#
#
###Stringers                     # C stringer dimentions
#HTh_str = Q_('0.025 m')            # height of the stringer
#HTw_str = Q_('0.025 m')            #width of the stringer
#HTt_str = Q_('0.003 m')            #thickness of the stringer
#
#
# z = 0.6696288358320435.4017773014992261                        #spanwise posotion in meters
# z *= Q_('meter')
# c = 0.0.0.0.0 dimensionless.0 dimensionless.1.1

##Ratio of height with respect to chord, airfoil coordinates
airfoilcoordinates = np.genfromtxt("../Airfoil.dat")    #Load coordinates
numberofcoordinates = np.size(airfoilcoordinates,0)  #Count total number of coordinates
airfoilinterpolant = sp.interpolate.interp1d(
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,0],
    airfoilcoordinates[0:int(numberofcoordinates/2)+1,1],kind = 'cubic') #Interpolate

perim_coords = np.genfromtxt("get_xy_from_perim.dat", skip_header=1)
perims = perim_coords[:,0]
x_cs = perim_coords[:,1]
y_cs = perim_coords[:,2]
perim_interpolant_x = sp.interpolate.interp1d(perims, x_cs, kind='cubic')
perim_interpolant_y = sp.interpolate.interp1d(perims, y_cs, kind='cubic')

perim_coords = np.genfromtxt("get_xy_from_perim_inv.dat", skip_header=1)
perims = perim_coords[:,0]
x_cs = perim_coords[:,1]
y_cs = perim_coords[:,2]
perim_interpolant_inv_x = sp.interpolate.interp1d(perims, x_cs, kind='cubic')
perim_interpolant_inv_y = sp.interpolate.interp1d(perims, y_cs, kind='cubic')
#print("xc=", perim_interpolant_inv_x(0.31415926535897932384626), perim_interpolant_inv_y(0.31415926535897932384626))

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

TR = airfoilordinate(0.15)*ChordR                            #max thickness root in m
TT = TR*t                                   #max thickness tip in m

def Chord_loc_Spar(zs,SparR,SparT):             #input spanwise location in m and
    ChSpar = SparR + (SparT-SparR)*(zs/s)  #Chord position of spar 1 with respect to leading edge
    return ChSpar

ChSpar1 = Chord_loc_Spar(z, Spar1R, Spar1T)
ChSpar2 = Chord_loc_Spar(z, Spar2R, Spar2T)

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
    angle = m.tan(dy/dx)
    angle *= Q_('rad')
    return angle

## Calculation of centriod

#Area of Spars
AreaSpar1 = HSpar1 * ThSpar1
AreaSpar2 = HSpar2 * ThSpar2

#Area additional clamps on spar 1
AreaClamps = ClampH * ClampW * 2                #Because there is a clamp on top as well as the bottom

## AREA OF STIFFENERS:
def calc_stringer_Area(w_str, h_str, t_str):
    b_1 = w_str / 2 - 0.5 * t_str
    h_1 = t_str

    b_2 = t_str
    h_2 = h_str

    b_3 = b_1
    h_3 = h_1

    # Calculating Area's of stringer components
    A_1 = b_1 * h_1
    A_2 = b_2 * h_2
    A_3 = A_1

    return (A_1 + A_2 + A_3)

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

A_stringer = calc_stringer_Area(w_str, h_str, t_str)                           #Area of one stringer
AreaStringers = A_stringer * N_stringers                #total area of stringer

## Area multiplied with the distance from the reference point (leading edge c=0)

#For the spars
AreaSpar1xc = AreaSpar1 * ChSpar1 * Chordlength
AreaSpar2xc = AreaSpar2 * ChSpar2 * Chordlength

#For the clamps
AreaClampsxc = AreaClamps * ChSpar1 * Chordlength

def get_xy_from_perim(perim_val, start_x=0, reverse=False):
    """
        NOTE: Special function for Tobias
        This function returns the x and y coordinates for a given perimeter


        Arguments: -perim_val (x coordinate)
                   -start_x: x_value from where the perimeter integration needs to start (IN RELATIVE COORD TO CHORD)
                   -dat_file_name (name of the airfoil data file)

        Returns: perim (Perimiter value)
    """

    #Air_data = np.genfromtxt(dat_file_name)  # Imports datapoints from airfoil data file
#
    #x_coords = Air_data[:81, 0]  # We only care about 1 half of the airfoil
    #x_coords = np.flip(x_coords, 0)  # Flip them so they are in a good order
    #y_coords = Air_data[:81, 1]  # We only care about 1 half of the airfoilfs
    #y_coords = np.flip(y_coords, 0)  # Flip them so they are in a good order
    #p = interp1d(x_coords, y_coords, kind ='cubic')  # Generate a poly spline based on the airfoil points

    perim = 0  # Set initial perimiter size to 0
    step = 0.001  # Step size for algorithm: increase will lead to faster computing times
    min_val = 10
    x_cs = np.arange(start_x, 1, step)
    if reverse == False:
        for i in range(len(x_cs)):
            perim += m.sqrt((step) ** 2 + (airfoilordinate(x_cs[i] + step) - airfoilordinate(x_cs[i])) ** 2)
            if abs(perim - perim_val) < min_val:
                x_coord = 0.5 * (x_cs[i] + x_cs[i] + step)
                y_coord = 0.5 * (airfoilordinate(x_cs[i]) + airfoilordinate(x_cs[i] + step))
                min_val = abs(perim - perim_val)
            else:
                break

        return(x_coord, y_coord)
    else:
        step *= -1
        for x_c in np.arange(start_x, 0, step):
            perim += m.sqrt((step) ** 2 + (airfoilordinate(x_c + step) - airfoilordinate(x_c)) ** 2)
            if abs(perim - perim_val) < min_val:
                x_coord = 0.5 * (x_c + x_c + step)
                y_coord = 0.5 * (airfoilordinate(x_c) + airfoilordinate(x_c + step))
                min_val = abs(perim - perim_val)
            else:
                break

        return (x_coord, -y_coord)

#tot_coords = np.array([0,0,0])
#for perim in np.arange(0,1,0.001):
#    x_coord, y_coord = get_xy_from_perim(perim, start_x=1, reverse=True)
#    coords = np.array([perim, x_coord, y_coord])
#    print(coords)
#    tot_coords = np.vstack((tot_coords, coords))
#np.savetxt("get_xy_from_perim_inv.dat", tot_coords, header="# XY coordinates calculated for a certain perimeters [dimensionless perim, dimensionless x_c, dimesnionless y_c]")

def get_perim_from_x(x_coor, inverse=False, dat_file_name="../Airfoil.dat"):
    """
    This function returns the perimeter value from the LE until the specified x-coordinate

    Arguments: x_coor (x coordinate)
               dat_file_name (name of the airfoil data file)

    Returns: perim (Perimiter value)
    """
    Air_data = np.genfromtxt(dat_file_name)  # Imports datapoints from airfoil data file

    x_coords = Air_data[:81, 0]  # We only care about 1 half of the airfoil
    x_coords = np.flip(x_coords, 0)  # Flip them so they are in a good order
    y_coords = Air_data[:81, 1]  # We only care about 1 half of the airfoil
    y_coords = np.flip(y_coords, 0)  # Flip them so they are in a good order
    p = interp1d(x_coords, y_coords, kind='cubic')  # Generate a poly spline based on the airfoil points

    perim = 0  # Set initial perimiter size to 0
    step = 0.0001  # Step size for algorithm: increase will lead to faster computing times
    # but lower accuracy
    # This for loop calculates the perimiter until the specified x-coordinate
    if inverse == False:
        for x_c in np.arange(0.0, x_coor, step):
            perim += np.sqrt((step) ** 2 + (p(x_c + step) - p(x_c)) ** 2)
    else:
        step *= -1
        for x_c in np.arange(1, x_coor, step):
            perim += np.sqrt((step) ** 2 + (p(x_c + step) - p(x_c)) ** 2)

    return perim

def lookup_xy_from_perim(norm_interp_x, norm_interp_y, inv_interp_x, inv_interp_y, perim, start_x_perim=0, inverse=False):
    if inverse == False:
        x_coor = norm_interp_x(perim+start_x_perim)
        y_coor = norm_interp_y(perim+start_x_perim)
        return (x_coor, y_coor)

    else:
        x_coor = inv_interp_x(perim+start_x_perim)
        y_coor = inv_interp_y(perim+start_x_perim)
        return (x_coor, y_coor)

# start_x_perim = get_perim_from_x(0.9, inverse=True)
# x_cs = np.array([])
# y_cs = np.array([])
# for per in np.arange(0.1, 0.89, 0.001):
#     x_coord, y_coord = lookup_xy_from_perim(perim_interpolant_x, perim_interpolant_y, perim_interpolant_inv_x, perim_interpolant_inv_y, per, start_x_perim=start_x_perim, inverse=True)
#     x_cs = np.append(x_cs, x_coord)
#     y_cs = np.append(y_cs, y_coord)
#
# plt.plot(x_cs, y_cs)
# plt.show()

def get_coord_from_perim(n_st, start_x, end_x, chord_l, dat_file_name="../Airfoil.dat"):
    """
    This function returns list of coordinate values where a stiffener is placed
    based on the spar locations and number of stiffeners. The stiffeners will
    be distributed equally along the perimiter section between the spars.

    Arguments: -n_st: number of stiffeners
               -start_x: x location of the LE spar
               -end_X: x location of the TE spar
               -dat_file_name: name of the airfoil data file

    Returns: List of x and y coordinates where the stiffeners are placed
    """
    Air_data = np.genfromtxt(dat_file_name)


    x_coords = Air_data[:81, 0]
    x_coords = np.flip(x_coords, 0)
    y_coords = Air_data[:81, 1]
    y_coords = np.flip(y_coords, 0)

    final_x_coords = np.array([])
    final_y_coords = np.array([])
    final_angles = np.array([])
    p = interp1d(x_coords, y_coords, kind='cubic')

    final_x_y_angle_coord = np.array([])
    perim = 0
    stif_perim = get_perim_from_x(end_x) - get_perim_from_x(start_x)
    perim_spacing = stif_perim / (n_st + 1)
    step = 0.0001
    i = 0
    for x_c in np.arange(start_x, end_x, step):
        perim += np.sqrt((step) ** 2 + (p(x_c + step) - p(x_c)) ** 2)
        if i >= n_st:
            break
        if abs(perim - perim_spacing) < 10e-5:
            x_coord = 0.5 * (x_c + x_c + step)
            y_coord = 0.5 * (p(x_c) + p(x_c + step))
            slope_angle = -(np.pi - np.arctan((p(x_c + step) - p(x_c)) / (step)))
            final_x_coords = np.append(final_x_coords, x_coord)
            final_y_coords = np.append(final_y_coords, y_coord)
            final_angles = np.append(final_angles, slope_angle)
            perim = 0
            i += 1

    # print(perim)
    #a_ran = np.arange(0, 1, 0.0001)
    #plt.plot([start_x, start_x], [0, p(start_x)], 'b')
    #plt.plot([end_x, end_x], [0, p(end_x)], 'b')
    #plt.plot(x_coords, y_coords, '+')
    #plt.plot(a_ran, p(a_ran), 'r')
    #plt.plot(final_x_coords, final_y_coords, 'go')
    #plt.axis((0, 1, 0, 1))
    #plt.show()
    return (final_x_coords*chord_l, final_y_coords*chord_l, final_angles), perim_spacing*chord_l

def stiffeners_centroid(x_y_angle_coords, h_str, w_str, t_str):

    x_coords = x_y_angle_coords[0]
    y_coords = x_y_angle_coords[1]
    angles = x_y_angle_coords[2]
    AX_cen = 0
    AY_cen = 0
    for i in range(len(x_coords)):
        x_cen = x_coords[i] + np.sin(angles[i])*h_str/2
        y_cen = y_coords[i] + np.cos(angles[i])*h_str/2
        AX_cen += x_cen*A_stringer
        AY_cen += y_cen*A_stringer
    X_cen = AX_cen/AreaStringers

    return(X_cen, 2*AX_cen)

def Area_Skin_x_c(Spar1, Spar2): #Input deminsionless chordwise location of spar 1 and spar 2
    n = 100 #number of sections
    dx = ((Spar2-Spar1)/n)
    x = Spar1
    Areaxc = Q_('0 m**3')
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

def length_Skin_x_c(Spar1, Spar2):                            #Input deminsionless chordwise location of spar 1 and spar 2
    n = 100 #number of sections
    dx = ((Spar2-Spar1)/n)
    x = Spar1
    arclength = 0
    for _ in range(n):
        x = x + dx
        dxlength = dx * Chordlength
        dylength = abs(airfoilordinate(x - dx) - airfoilordinate(x)) * Chordlength
        dlength = np.sqrt(dxlength**2+dylength**2)
        arclength += dlength
    return arclength

## Area of cell enclosed by the wing skin and the stringers
def Area_cell(ChSpar2, ChSpar1):
    n = 100 #number of sections
    dx = ((ChSpar2-ChSpar1)/n)
    area_cell = 0
    x = ChSpar1
    for i in range(n):
        x = x + dx
        dxlength = dx * Chordlength
        area_cell = dxlength*airfoilordinate(x)*Chordlength
        area_cell = area_cell*2
    return area_cell

area_cell = Area_cell(ChSpar2, ChSpar1)
skin_length = length_Skin_x_c(ChSpar1, ChSpar2) # skin length of top part of airfoil (dimensionless)

x_y_angle_coords, perim_spacing = get_coord_from_perim(N_stringers/2, ChSpar1, ChSpar2, Chordlength)
X_cen_strs, A_stringer_x_c = stiffeners_centroid(x_y_angle_coords, h_str, w_str, t_str)
Area_x_c = AreaSpar1xc + AreaSpar2xc + Area_Skin_x_c(ChSpar1, ChSpar2) + A_stringer_x_c + AreaClampsxc
Area = AreaSpar1 + AreaSpar2 + Area_Skin(ChSpar1, ChSpar2)+ AreaStringers + AreaClamps

centroidspars= (AreaSpar1xc + AreaSpar2xc)/(AreaSpar1 + AreaSpar2)
centroidskin=Area_Skin_x_c(ChSpar1, ChSpar2)/Area_Skin(ChSpar1, ChSpar2)
centroidstringer= A_stringer_x_c/AreaStringers

# print('spars', centroidspars/Chordlength)
# print('skin', centroidskin/Chordlength)
# print('stringer', centroidstringer/Chordlength)
centroid = (Area_x_c/Area)/Chordlength
centroidlength = Area_x_c/Area

print(centroid)
# print(ChSpar1)

#
# ''''Calculate volume of the '''
#
# def Vol_wingbox(Spar1, Spar2, Chordlength):                            #Input deminsionless chordwise location of spar 1 and spar 2
#     n = 100 #number of sections
#     dx = ((Spar2-Spar1)/n)
#     x = Spar1
#     Area = 0
#     for i in range(n):
#         y = airfoilordinate(x)
#         Area = Area + y*dx*(Chordlength**2)
#         x = x + dx
#     Area = Area * 2                                 # Area of both sides of the airfoil
#     return Area
#
# nx = 30
# tankstop = Q_('1 m')
# z = Q_('0 m')
# Voltank = Q_('0 m**3')
#
# while z < tankstop:
#     Areatank = Vol_wingbox(Chord_loc_Spar(z, Spar1R, Spar1T), Chord_loc_Spar(z, Spar2R, Spar2T), length_chord(z))
#     section = tankstop/nx
#     Voltank = Voltank + Areatank * section
#     z = z + section
# print('Voltank', Voltank)
#
print('z', z)