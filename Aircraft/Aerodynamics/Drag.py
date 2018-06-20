# -*- coding: utf-8 -*-
"""
Name: DATCOM Drag
Department: Aerodynamics
Last updated: 13/06/2018 11:00 by Emma
"""

#Imports
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
import math as m
import numpy as np
import scipy as sp
from Geometry import Geometry
from pint import UnitRegistry
import Aeroprops



#Input variables, to be filled in by user


AoA = Q_("2 deg") #fill in angle of attack flight condition drag to be known
AoA = AoA.to(ureg.rad)
Rwb = 1.05#Read in DATCOM p1164, wing body interference factor
Rhtb = 1.05
Rvtb= 1.05
CL_trim = 0#lift coefficient needed to trim aircraft

height = Q_('500 m') # height in m
velocity = Q_('30  m/s') #velocity in m/s
loc_max_tc_wing = 0 # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
loc_max_tc_ht = 0 # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
loc_max_tc_vt = 0 # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
Sref = Geometry.Wing.S
CL_wing = Aeroprops.CL_alpha_wing * AoA #max lift coefficient
AR_wing = Geometry.Wing.A
ih = Geometry.H_tail.i_h #incidende angle HT
CL_ht = Aeroprops.CL_alpha_ht * AoA.to(ureg.rad) - Aeroprops.de_da * \
                AoA.to(ureg.rad) + ih#max lift coefficient ht
AR_ht = Geometry.H_tail.A
span_eff = 0.792 #span efficiency factor for lift induced drag calcs
length_fus = Geometry.Fuselage.l_f #fuselage length
fus_cs_area = Geometry.Fuselage.front_A
S_wet_fus = Geometry.Fuselage.S_wet_f
Sht = Geometry.H_tail.S
Svt = Geometry.V_tail.S



#Standard values Inputs


density_0 = Q_('1.225 kg/m**3')
Temp0 = Q_('288.15 K')
g = Q_('9.80665 m / s**2')
R = Q_('287.05 m**2 * K / s**2')
lam = Q_('0.0065 K / m')
Temp = Temp0 - lam * height
Rls_para = 1.07 #obtained from DATCOM, won't change with slight aircraft redesign


#Input variables, linked to aircraft geometry files etc.


chord_wing = Geometry.Wing.c_avg #average chord
chord_ht = Geometry.H_tail.MAC
chord_vt = Geometry.V_tail.MAC
airfoil_co_wing = np.genfromtxt("../airfoil.dat")
airfoil_co_ht = np.genfromtxt("../airfoil.dat")
airfoil_co_vt = np.genfromtxt("../airfoil.dat")


#general calculations outside of for loop
viscosity = 1.458e-6 * Temp.magnitude**1.5 * (1 / (Temp.magnitude + 110.4))
viscosity *= Q_('1 N * s / m**2')
density = density_0 * (Temp / Temp0)**((g.magnitude/ lam.magnitude / R.magnitude) - 1) #ISA formula


#Interpolation of graphs
#cut off Reynolds number
roughness_points = np.log10(np.array([(1.9e3,1e5),(6e5,1e4),(1e6,1.8e4),
                                          (7e6,9e4),(1e8,1.2e6),(9e8,1e7)]))
cutoffre = roughness_points[:,0] #x values selection
adm_rough = roughness_points[:,1] # y values selection, logarithmic scale
    
#skin friction coefficient
Cf_points = np.array([(3.5e5,0.0055),(5e5,0.005),(1.7e6,0.004),(9.5e6,0.003),(1.5e8,0.002),
                          (1e9,0.0016)])
inputre = np.log10(Cf_points[:,0]) #xvals, logaritmhic scale
outputcf = Cf_points[:,1] #yvals
CD0_array = np.array([])   
CDi_array = np.array([])

for lift_surface in np.array(['wing', 'ht', 'vt']):
    if lift_surface == 'wing':
        chord = chord_wing
        loc_max_tc = loc_max_tc_wing
        airfoil_co = airfoil_co_wing
        Swet = 2 * Sref
        CL = CL_wing
        AR = AR_wing
        e = span_eff
    elif lift_surface == 'ht':
        chord = chord_ht
        loc_max_tc = loc_max_tc_ht
        airfoil_co_ht = airfoil_co_ht
        Swet = 2 * Sht
        CL = CL_ht
        AR = AR_ht
        e = 0.5
    elif lift_surface == 'vt':
        chord = chord_vt
        loc_max_tc = loc_max_tc_vt
        airfoil_co_vt = airfoil_co_vt
        Swet = 2 * Svt
        CL = 0
        e = 0.5 #value doesn't matter
        AR = 1 #value doesn't matter for this, as long as CL = 0
        
    #(Reynolds number calculation        
    Reynolds = density * velocity * chord.to(ureg.meter) / viscosity #definition Reynolds number
    print('Reynolds',Reynolds)
    # Interpolation for Cutoff Reynolds number
    chordinch = chord.to(ureg.inch) #convert chord to inch for graph implementationn logarithmic scale
    print('chordinch', chordinch)
    input_cut_re = np.log10(chordinch.magnitude / 0.3e-3)  #x-input for interpolation graph
    cutoff_Reynolds = 10**sp.interp(input_cut_re, adm_rough, cutoffre) #linear fit
    print('cutoffRE',cutoff_Reynolds)
    
    Re_used = min(Reynolds, cutoff_Reynolds) #use lowest Reynolds number
    print('reusde', Re_used)
    
    # Skin friction coefficient interpolation and calculation
    C_skinfric = (sp.interp(np.log10(Re_used), inputre,outputcf))
    print('skinfricco',C_skinfric)
    
    #calculate average t/c airfoil
    upper_half = airfoil_co[0:int((len(airfoil_co)+1)/2),1]
    average_thickness = 2 * np.average(upper_half[:]) #*2 since symmetric
    print('thick',average_thickness)
    
    #L parameter DATCOM
    if loc_max_tc == 0:
        Lpara = 2
    elif loc_max_tc == 1:
        Lpara = 1.2
        
        
    #cdo calculations
    CDO_lifting_s = C_skinfric * (1 + Lpara * average_thickness + 100 *\
                                  average_thickness ** 4) * Rls_para * \
                                  Swet/Sref
    CD0_array = np.append(CD0_array, CDO_lifting_s )
    print('CDO of', lift_surface, '=', CDO_lifting_s)
    #Lift induced drag
    CDi = CL**2 / m.pi / AR / span_eff
    CDi_array = np.append(CDi_array, CDi)
    print('CDi of', lift_surface,'=', CDi)


# Body Drag

#Reynolds number calculation        
Reynolds_fus = density * velocity * length_fus / viscosity #definition Reynolds number

# Interpolation for Cutoff Reynolds number

length_fus_inch = length_fus.to(ureg.inch) #convert chord to inch for graph implementationn logarithmic scale
input_cut_re_fus = np.log10(length_fus_inch.magnitude / 0.25e-3)  #x-input for interpolation graph
cutoff_Reynolds_fus = 10**sp.interp(input_cut_re_fus, adm_rough, cutoffre) #linear fit

Re_used_fus = min(Reynolds_fus, cutoff_Reynolds_fus) #use lowest Reynolds number
print('reusedfus', Re_used_fus)

# Skin friction coefficient interpolation and calculation
C_skinfric_fus = sp.interp(np.log10(Re_used_fus), inputre,outputcf)
print('skinfric fus', C_skinfric_fus)

# equivalent fuselage diameter to us rotational body calculation DATCOM
dia_eq_fus = np.sqrt(fus_cs_area / 0.7854)
print(dia_eq_fus, 'dia eq fus')
fin_ratio = length_fus / dia_eq_fus #fuselage fineness ratio
print('fineness',fin_ratio)

#zero lift body drag

CD0_fus = C_skinfric_fus * (1  +  60 / (fin_ratio**3)  +  0.0025 * fin_ratio)  * \
S_wet_fus / Sref

print('fuselage zero lift drag', CD0_fus)


#Total zero lift drag
CD0_tot = CD0_array [0] * Rwb + CD0_fus * Rwb * Rhtb * Rvtb + CD0_array[1] *\
 Rhtb  + CD0_array[2] * Rvtb
print('cdotot', CD0_tot)


#trim drag
downwash = Aeroprops.de_da * AoA
Delta_CD_trim = ((CD0_array[1] + CL_trim**2 / m.pi / AR_ht / 0.5) * \
                 np.cos(downwash) + CL_trim * np.sin(downwash)) * \
                 Sht / Sref * Aeroprops.q_qinf_ratio
                 
        
#miscellaneous drag
CD_canopy = Q_('0.04 m**-2') * Geometry.Fuselage.A_max_canopy
print(CD_canopy,'cap')
CD_lg = Q_('0.458 m**-2') * Geometry.Landing_gear.lg_wheel_d * Geometry.Landing_gear.lg_wheel_w
print('landing gear',CD_lg)


#total drag calculation

CD_tot = CD0_tot + CD_canopy + CD_lg + Delta_CD_trim + CDi_array[0] + CDi_array[1]  

Drag = 0.5 * density * velocity**2 * CD_tot * Sref
    

#printing stuff to reference file
with open('Aeroprops.py', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

data[20] = 'CD0_wing = Q_(\"' + str(CD0_array[0]) + '\")\n'
data[21] = 'CD0_ht = Q_(\"' + str(CD0_array[1]) + '\")\n'
data[22] = 'CD0_vt = Q_(\"' + str(CD0_array[2]) + '\")\n'
data[23] = 'CD_canopy = Q_(\"' + str(CD_canopy) + '\")\n'
data[24] = 'CD_lg = Q_(\"' + str(CD_lg) + '\")\n'
data[25] = 'CD0_tot = Q_(\"' + str(CD0_tot) + '\")\n'

# and write everything back
with open('Aeroprops.py', 'w') as file:
    file.writelines(data)







