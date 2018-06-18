# -*- coding: utf-8 -*-
"""
Name: Slat Force
Department: Aerodynamics
Last updated: 13/06/2018 14:14 by Emma
"""

#imports
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
import numpy as np
import math as m
from Geometry import Geometry

#Manual inputs
#chordlength = .05 #determine part of chord used for slats
slatwidth = Geometry.Wing.b  #width of the slats
h = Q_('100 m') #altitude of flight
Velocity = Q_(' 30 m/s') #aircraft velocity at which slats are deployed
GMAC = Geometry.Wing.c_avg
mass_sys = 20 #mass of slat system
#%% Slat sizing for optimal max lift increase

# optimization arrays:
deflection = np.linspace(15,25,20) #deflection in degrees
deflection_rad = np.radians(deflection)
slatchordratio = np.linspace(0,0.15,20) #dimensionless
verticaldeflection = np.linspace(0.05,0.0845,20) #in meter

#DATCOM graph interpolation
maxLeffect_datapoints = np.array([(0,0),(.05,.84),(.1,1.2),(.15,1.31),(.2,1.6),\
                                  (.25,1.75),(.3,1.825)])
maxL_x = maxLeffect_datapoints[:,0]
maxL_y = maxLeffect_datapoints[:,1]
eta_max = .8
etadelta_datapoints = np.array([(15,1),(17.5,.95),(20,.9),(22.5,.82),(25,7.6),(30,.59)])
ed_x = etadelta_datapoints[:,0]
ed_y = etadelta_datapoints[:,1]

#optimization loop
optimum = 0
darray = np.array([])
carray = np.array([])
garray = np.array([])
for defl in deflection_rad:
    for c_slat in slatchordratio:
        for gap in verticaldeflection:
            maxlifteffect = np.interp(c_slat,maxL_x,maxL_y)
            da_corr_fact = np.interp(defl,ed_x,ed_y)
            delta_clmax = maxlifteffect * eta_max * da_corr_fact * defl *\
            (GMAC.magnitude + gap)/GMAC.magnitude
            
            if delta_clmax > optimum:
                optimum = delta_clmax
                darray = np.append(darray,defl)
                carray = np.append(carray,c_slat)
                garray = np.append(garray,gap)


defl_angle = darray[-1]
chordslat = carray[-1]
gapslat = garray[-1]
print('best delta_clmax =', optimum)
print('deflection angle optimum =', m.degrees(defl_angle) )
print('chord optimum = ', chordslat)
print('gap optimum = ', gapslat)              

#%% Calculate surface manin airfoil that will be slat (length in m)
airfoildata = np.genfromtxt("../airfoil.dat") #import airfoil
halfairfoil = airfoildata[0:81,:] # take upper half of airfoil coordinates
counter = -1
distance = np.array([])
for xco in halfairfoil[:,0]:
    counter = counter + 1
    if xco <= chordslat:
        #print('xco', halfairfoil[counter,0])
        new_d = m.hypot(halfairfoil[counter,0] - halfairfoil[counter - 1,0],\
                        halfairfoil[counter,1] - halfairfoil[counter - 1,1])
        distance = np.append(distance,new_d)
        #print(counter)
print('distance',distance)
slatlength = np.sum(distance) * Geometry.Wing.c_avg #upper part of slat airfoil length
print('slatlength',slatlength)


#ISA calculations
density_0 = Q_('1.225 kg/m**3')
Temp0 = Q_('288.15 K')
g = Q_('9.80665 m / s**2')
R = Q_('287.05 m**2 * K / s**2')
lam = Q_('0.0065 K / m')
p0 = Q_('101325 N/m**2')
Temp = Temp0 - lam * h
density = density_0 * (Temp / Temp0)**((g.magnitude/ lam.magnitude / R.magnitude) - 1)


# Slat force calculation
slatarea = slatlength * slatwidth
print('area',slatarea)
P = 0.5 * density * Velocity**2 * slatarea #F = pressure times area
print('Slat force required =', P)

#suction peak calculation
cp_file = np.genfromtxt('../Cp_slat.dat')
cp_array = np.array([])
for coor in np.arange(0,len(cp_file)):
    if cp_file[coor,0] <= chordslat and cp_file[coor,1] > 0:
        cp_array = np.append(cp_array, cp_file[coor,2])
    
cp_ave = np.average(cp_array)
F = 0.5 * density * Velocity**2 * cp_ave * slatarea

acc = F.magnitude / mass_sys
time = m.sqrt(gapslat/-acc)


#%% electric actuator size + stroke length
stroke = gapslat / m.cos(defl_angle)
print('stroke length',stroke)


total_length = 0.2
rod = np.linspace(0.07,.08,50)
pushrod = stroke + 0.02
rodangle = np.radians(np.linspace(0,90,50))
for r in rod:
    h = r - pushrod * m.sin(defl_angle)
    l = total_length - m.cos(defl_angle) * pushrod
    midrod = m.sqrt(h**2 + l**2)
    for a in rodangle:
        new_h = r * m.cos(a)
        new_x = r * m.sin(a)
        midrodh = new_h - (pushrod - stroke) * m.sin(defl_angle)
        midrodx = m.sqrt(midrod**2 - midrodh**2)
        gain = midrodx - l + new_x
        if gain > gapslat:
            print('for a gain of',gain, 'a rod angle of', m.degrees(a), ' and a rod \
                  length of ', r, ' is need')
            
        

#%% Forces in the fuselage

#principle of work to do forces in fuselage

rodangle = m.radians(86.3265)#choose from a printed value above
work_slat = P * stroke / 2 #work = force times displacement 2 actuator for this force
work_slat *= Q_('N * m')
F_act = Q_('500 N') #electric actuator force
displ_fus = work_slat / F_act
print('linear displacement rod =', displ_fus)
fus_rod_length = (2 * (displ_fus * m.sin(rodangle))**2)
print('rod in fuselage has a length of', fus_rod_length)
rod2 = fus_rod_length - fus_rod_length * m.cos(rodangle/2)
print('small rod = ', rod2)





"""
x_length = (chordlength - 0.01) * Geometry.Wing.c_avg.magnitude
print('x',x_length)
y_length = (Geometry.Wing.T_Cmax / 2 - 0.03) * Geometry.Wing.c_avg.magnitude
print('y',y_length)
angle = m.degrees(m.atan(y_length/x_length))
print(angle,' = angle')
stroke = m.sqrt(x_length**2 + y_length**2)
print('stroke length =', stroke)

x_attach = (0.184 - 0.2 * chordlength ) * Geometry.Wing.c_avg.magnitude
y_attach = 2/3 * Geometry.Wing.T_Cmax / 2 * Geometry.Wing.c_avg.magnitude
angle_actuator = m.degrees(m.atan(y_attach / x_attach))
print('actuator anagl', angle_actuator)
actuator = m.sqrt(x_attach**2 + y_attach**2)
print('actuator length',actuator)
"""
