# -*- coding: utf-8 -*-
"""
Name: DATCOM Drag
Department: Aerodynamics
Last updated: 08/06/2018 16:53 by Emma
"""

#Imports

import math as m
import numpy as np
import scipy as sp
from Geometry import Geometry
from pint import UnitRegistry
ureg = UnitRegistry()



#Input variables, to be filled in by user


Temp = 100# Temperature in Kelvin
Temp *= Q_('K')
velocity = Q_('10  m/s') #velocity in m/s
loc_max_tc_wing = # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
loc_max_tc_ht = # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
loc_max_tc_vt = # 0 if t/c max is < 30%c, 1 if t/c max is >= 30%c
Sref = Geometry.Wing.S
CL_wing = 
AR_wing = Geometry.Wing.A
CL_ht =
AR_ht = Geometry.H_tail.A
span_eff = 0.792 #span efficiency factor for lift induced drag calcs


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
airfoil_co_ht =
airfoil_co_vt = 


#general calculations outside of for loop
viscosity = 1.458e-6 * Temp**1.5 * (1 / (Temp + 110.4))
viscosity *= Q_('N * s / m**2')
density = density_0 * (Temp / Temp0)**((g/ lam / R) - 1) #ISA formula




for lift_surface in np.array(['wing', 'ht', 'vt']):
    if lift_surface == 'wing':
        chord = chord_wing
        loc_max_tc = loc_max_tc_wing
        airfoil_co = airfoil_co_wing
        Swet = 2 * Geometry.Wing.A
        CL = CL_wing
        A = AR_wing
    elif lift_surface == 'ht':
        chord = chord_ht
        loc_max_tc = loc_max_tc_ht
        airfoil_co_ht = airfoil_co_ht
        Swet = 2 * Geometry.H_tail.A
        CL = CL_ht
        A = AR_ht
    elif lift_surface == 'vt':
        chord = chord_vt
        loc_max_tc = loc_max_tc_vt
        airfoil_co_vt = airfoil_co_vt
        Swet = 2 * Geometry.V_tail.A
        CL = 0
        A = 1 #value doesn't matter for this, as long as CL = 0
        
    #(Reynolds number calculation        
    Reynolds = density * velocity * chord / viscosity #definition Reynolds number
    
    # Interpolation for Cutoff Reynolds number
    roughness_points = np.log10(np.array([(1.9e3,1e5),(6e5,1e4),(1e6,1.8e4),
                                          (7e6,9e4),(1e8,1.2e6),(9e8,1e7)]))
    cutoffre = roughness_points[:,0] #x values selection
    adm_rough = roughness_points[:,1] # y values selection, logarithmic scale
    chordinch = chord.ito(ureg.inch) #convert chord to inch for graph implementationn logarithmic scale
    input_cut_re = np.log10(chordinch / 0.3e-3)  #x-input for interpolation graph
    cutoff_Reynolds = 10**sp.interp(input_cut_re, adm_rough, cutoffre) #linear fit
    
    Re_used = min(Reynolds, cutoff_Reynolds) #use lowest Reynolds number
    
    
    # Skin friction coefficient interpolation and calculation
    Cf_points = np.array([(5e5,0.005),(1.7e6,0.004),(9.5e6,0.003),(1.5e8,0.002),
                          (1e9,0.0016),(3.5e5,0.0055)])
    inputre = np.log10(Cf_points[:,0]) #xvals, logaritmhic scale
    outputcf = Cf_points[:,1] #yvals
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
    
    #Lift induced drag
    CDi = CL**2 / m.pi / AR / span_eff

# Body Drag



    
    








