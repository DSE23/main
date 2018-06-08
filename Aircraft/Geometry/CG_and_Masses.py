"""
Created on Thu Jun  7 14:34:37 2018

@author: jurian

StefX CG locations + masses
"""
import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

W_wing = Q_("121 kg")                # [kg] Mass of the wing
W_htail = Q_("18 kg")                # [kg] Mass of H_tail
W_vtail = Q_("6 kg")                 # [kg] Mass of V_tail
W_fus = Q_("82 kg")                  # [kg] Mass of Fuselage
W_gear = Q_("58 kg")                 # [kg] Mass of landing gear
W_engine = Q_("324 kg")              # [kg] Mass of engine
W_prop = Q_("0 kg")                  # [kg] Mass of propellor
W_fuelsys = Q_("10 kg")              # [kg] Mass of fuel system
W_hydraulic = Q_("1 kg")             # [kg] Mass of hydraulics
W_flightcontrol = Q_("20 kg")        # [kg] Mass of flight control
W_avionics = Q_("17 kg")             # [kg] Mass of Avionics
W_elecsys = 