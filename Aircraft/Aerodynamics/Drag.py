# -*- coding: utf-8 -*-
"""
Name: DATCOM Drag
Department: Aerodynamics
Last updated: 08/06/2018 16:53 by Emma
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
ureg = UnitRegistry()
from Aeroprops import Aeroprops



#Input variables, to be filled in by user


AoA = 1#fill in angle of attack flight condition drag to be known
Rwb = 1#Read in DATCOM p1164, wing body interference factor

Temp = 100# Temperature in Kelvin
Temp *= Q_('K')
velocity = Q_('10  m/s') #velocity in m/s
loc_max_tc_wing = 0 # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
loc_max_tc_ht = 0 # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
loc_max_tc_vt = 0 # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
Sref = Geometry.Wing.S
CL_wing = Aeroprops.CL_alpha_wing * AoA.ito(ureg.rad) #max lift coefficient
AR_wing = Geometry.Wing.A
ih = Geometry.H_tail.i_h #incidende angle HT
CL_ht = Aeroprops.CL_alpha_ht * AoA.ito(ureg.rad) - Aeroprops.de_da * 
                AoA.ito(ureg.rad) + ih#max lift coefficient ht
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
Rls_para = 1.07 #obtained from DATCOM, won't change with slight aircraft redesign


#Input variables, linked to aircraft geometry files etc.


chord_wing = Geometry.Wing.c_avg #average chord
chord_ht = Geometry.H_tail.MAC
chord_vt = Geometry.V_tail.MAC
airfoil_co_wing = np.genfromtxt("../airfoil.dat")
airfoil_co_ht = np.genfromtxt("../airfoil.dat")
airfoil_co_vt = np.genfromtxt("../airfoil.dat")


#general calculations outside of for loop
viscosity = 1.458e-6 * Temp**1.5 * (1 / (Temp + 110.4))
viscosity *= Q_('N * s / m**2')
density = density_0 * (Temp / Temp0)**((g/ lam / R) - 1) #ISA formula


#Interpolation of graphs
#cut off Reynolds number
roughness_points = np.log10(np.array([(1.9e3,1e5),(6e5,1e4),(1e6,1.8e4),
                                          (7e6,9e4),(1e8,1.2e6),(9e8,1e7)]))
    cutoffre = roughness_points[:,0] #x values selection
    adm_rough = roughness_points[:,1] # y values selection, logarithmic scale
    
#skin friction coefficient
Cf_points = np.array([(5e5,0.005),(1.7e6,0.004),(9.5e6,0.003),(1.5e8,0.002),
                          (1e9,0.0016),(3.5e5,0.0055)])
    inputre = np.log10(Cf_points[:,0]) #xvals, logaritmhic scale
    outputcf = Cf_points[:,1] #yvals
    
for lift_surface in np.array(['wing', 'ht', 'vt']):
    if lift_surface == 'wing':
        chord = chord_wing
        loc_max_tc = loc_max_tc_wing
        airfoil_co = airfoil_co_wing
        Swet = 2 * Sref
        CL = CL_wing
        A = AR_wing
    elif lift_surface == 'ht':
        chord = chord_ht
        loc_max_tc = loc_max_tc_ht
        airfoil_co_ht = airfoil_co_ht
        Swet = 2 * Sht
        CL = CL_ht
        A = AR_ht
    elif lift_surface == 'vt':
        chord = chord_vt
        loc_max_tc = loc_max_tc_vt
        airfoil_co_vt = airfoil_co_vt
        Swet = 2 * Svt
        CL = 0
        A = 1 #value doesn't matter for this, as long as CL = 0
        
    #(Reynolds number calculation        
    Reynolds = density * velocity * chord / viscosity #definition Reynolds number
    
    # Interpolation for Cutoff Reynolds number
    chordinch = chord.ito(ureg.inch) #convert chord to inch for graph implementationn logarithmic scale
    input_cut_re = np.log10(chordinch / 0.3e-3)  #x-input for interpolation graph
    cutoff_Reynolds = 10**sp.interp(input_cut_re, adm_rough, cutoffre) #linear fit
    
    Re_used = min(Reynolds, cutoff_Reynolds) #use lowest Reynolds number
    
    
    # Skin friction coefficient interpolation and calculation
    C_skinfric = np.log10(sp.interp(np.log10(Re_used), inputre,outputcf))
    
    
    #calculate average t/c airfoil
    upper_half = airfoil_co[0:int((len(airfoil_co)+1)/2),1]
    average_thickness = 2 * np.average(upper_half[:]) #*2 since symmetric
    
    
    #L parameter DATCOM
    if loc_max_tc == 0:
        Lpara = 2
    elif loc_max_tc == 1:
        Lpara = 1.2
        
        
    #cdo calculations
    CDO_lifting_s = C_skinfric * (1 + Lpara * average_thickness + 100 *
                                  average_thickness ** 4) * Rls_para *
                                  Swet/Sref
    print('CDO of', lift_surface, '=', CDO_lifting_s)
    #Lift induced drag
    CDi = CL**2 / m.pi / AR / span_eff
    print('CDi of', lift_surface,'=', CDi)


# Body Drag

#Reynolds number calculation        
Reynolds_fus = density * velocity * length_fus / viscosity #definition Reynolds number

# Interpolation for Cutoff Reynolds number

length_fus_inch = length_fus.ito(ureg.inch) #convert chord to inch for graph implementationn logarithmic scale
input_cut_re_fus = np.log10(length_fus_inch / 0.3e-3)  #x-input for interpolation graph
cutoff_Reynolds_fus = 10**sp.interp(input_cut_re_fus, adm_rough, cutoffre) #linear fit

Re_used_fus = min(Reynolds_fus, cutoff_Reynolds_fus) #use lowest Reynolds number


# Skin friction coefficient interpolation and calculation
C_skinfric_fus = np.log10(sp.interp(np.log10(Re_used_fus), inputre,outputcf))

# equivalent fuselage diameter to us rotational body calculation DATCOM
dia_eq_fus = m.sqrt(fus_cs_area / 0.7854)
fin_ratio = length_fus / dia_eq_fus #fuselage fineness ratio

#zero lift body drag

CDO_fus = C_skinfric_fus * (1  +  60 / fin_ratio**3  +  0.0025 * fin_ratio)  *  
            S_wet_fus / fus_cs_area

print('fuselage zero lift drag', CDO_fus)


#Total zero lift drag
CDO_tot = CD0_wing * Rwb + CD0_fus * Rwb * Rhtb * Rvtb + CD0_ht * Rhtb  + CD0_vt * Rvtb


#trim drag
downwash = Aeroprops.de_da * AoA
Delta_CD_trim = ((CD0_ht + CDi_ht) * m.cos(downwash) + CL_ht * m.sin(downwash)) *
            Sht / Sref * Aeroprops.q_qinf_ratio
    
    








