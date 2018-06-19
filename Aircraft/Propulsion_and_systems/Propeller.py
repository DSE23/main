"""
Name: Propeller
Department: Propulsion and Aircraft Systems
Last updated: 06/06/2018 10:25 by Ties
"""

import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

import numpy as np
from Geometry import Geometry

# Returns the propeller efficiency, thrust and pitch at the hub for the maximum power value, for a given airspeed inbetween 15m/s and 200m/s.

def Thrustcalc(V0):
    #Required values
    D = Geometry.Prop.Diameter.magnitude          #Diameter of the propeller
    R = D/2                             #Radius of the propeller
    Rhub = 0.20                         #Radius of the hub. Optimal at 0.31
    Elements = 10                     #Number of elements on the blade
    RPM = 2700                          #constant speed prop
    Omega = RPM/60*2*np.pi              #rotational speed, rad/s
    TotalTwist = 33.6/180*np.pi         #Total twist of the propeller blade
    PitchRange = np.arange(35,70,0.1)   #Range of starting pitch of the hub which will be tested.
    rho = 1.225                         #Density at SL
    Nprop = 4                           #Number of propellers

    running = True
    Thrustlst = []
    Torquelist = []
    efflist = []
    BetaHublst = []
    Rvel = []

    for i in PitchRange:
        BetaHub = i / 180 * np.pi
        BetaTip = BetaHub - TotalTwist

        w = np.ones(Elements) * 0.001
        BetaHublst.append(BetaHub)

        # iterative starting values
        Alphai = 0

        # Dimensions
        deltar = (R - Rhub) / Elements                          # Element length
        r = Rhub + deltar * np.arange(Elements) + 0.5 * deltar  # Middle of elements in array
        x = r / R                                               # percentage of span
        c = 0.152 - 0.05*x                                    # Element cord, an estimation from some example. Hopefully MT will respond and we will get real values.

        deltaA = c * deltar                                     # Element area

        # Speeds
        Propspeed = Omega * r                                   # Velocity from rotation for every element
        VR = np.sqrt(V0 ** 2 + Propspeed ** 2)
        VE = VR * np.cos(Alphai)
        
        
        count = 0

        # flow angels
        phi = np.arctan(V0 / Propspeed)
        beta = BetaHub * (1 - (np.arange(Elements) + 0.5) / Elements) + BetaTip * (np.arange(Elements) + 0.5) / Elements
        AlphaZL = -3.80 * np.pi / 180

        Iterating = True

        while Iterating:
            Alphai = np.arctan(w / VE)
            Alpha = beta - phi - Alphai + AlphaZL
            AlphaDeg = Alpha * 180 / np.pi

            ClProp = 1E-07 * AlphaDeg ** 6 - 4E-06 * AlphaDeg ** 5 + 2E-05 * AlphaDeg ** 4 + 0.0003 * AlphaDeg ** 3 - 0.0019 * AlphaDeg ** 2 + 0.1119 * AlphaDeg + 0.456
            CdProp = -3E-09 * AlphaDeg ** 6 + 6E-08 * AlphaDeg ** 5 + 5E-06 * AlphaDeg ** 4 - 1E-04 * AlphaDeg ** 3 + 0.0005 * AlphaDeg ** 2 + 0.0011 * AlphaDeg + 0.0043
            dL = 0.5 * rho * ClProp * deltaA * VE ** 2
            dD = 0.5 * rho * CdProp * deltaA * VE ** 2

            dT = dL * np.cos(phi + Alphai) - dD * np.sin(phi + Alphai)
            dQ = r * (dL * np.sin(phi + Alphai) + dD * np.cos(phi + Alphai))
            dP = dQ * Omega

            Thrust = sum(dT) * Nprop
            Torque = sum(dQ) * Nprop
            Power = sum(dP) * Nprop

            CT = Thrust / (rho * (RPM / 60) ** 2 * D ** 4)
            CQ = Torque / (rho * (RPM / 60) ** 2 * D ** 5)
            CP = Power / (rho * (RPM / 60) ** 3 * D ** 5)

            J = V0 / (RPM / 60 * D)

            VE = np.sqrt((w + V0) ** 2 + Propspeed ** 2)
            Fw = 8 * np.pi * r / Nprop / c * w - VE / (w + V0) * (ClProp * Propspeed - CdProp * (w + V0))
            Fwaccent = 8 * np.pi * r / Nprop / c - ClProp * Propspeed * (1 / VE - VE / (V0 + w) ** 2) + CdProp * (w + V0) / VE

            wtemp = w - Fw / Fwaccent

            if max(abs(wtemp - w)) < 0.001:
                Iterating = False

            if count > 200:
                Iterating = False

            w = wtemp
            count = count + 1
        etaP = J * CT / CP*100

        if Power >= 235000 and Power <= 250000 and running:
            running = False
            Final = [efflist[-1],Thrustlst[-1],i,w,Alphai,VR[-Elements:-1]]
        efflist.append(etaP)
        Thrustlst.append(Thrust)
        Torquelist.append(Torque)
        Rvel.append(VR)
    return Final


#Returns propeller efficiency, Thrust, Pitchangle, induced velocity, and Induced AOA for max power.
# For Alphai: w cos(Alphai+phi) is in airflow direction. w sin(Alphai + phi) is naar buiten.


D = 1.90            #Diameter of the propeller
R = D/2             #Radius of the propeller
Rhub = 0.20         #Radius of the huboptimal at 0.31
Elements = 1000
P = 235000
rho = 1.225

Tstatic = 0.85*P**(2/3)*(2*rho*R**2*np.pi)**(1/3)*(1-Rhub**2/(R**2))    #Maximum static thrust.



# Coordinate system:
# Origin is in propeller attachment to engine
# X axis: parallel to the crankshaft centerline, pointing towards tail
# Y axis: up
# Z axis: left

# Start defining global variables for easy editing from elsewhere
# For explanations of the variables defined here, see below, where they are given values


def initialise_mass(inp):
    global mass
    mass = inp


def initialise_propxcg(inp):
    global xcg
    xcg = inp


def initialise_propinertia_y(inp):
    global iyg
    iyg = inp


def initialise_propinertia_y(inp):
    global iyg
    iyg = inp


# End defining global variables

# Start assigning values to variables
# Propeller mass

# CG


# Mass moments of inertia
ixg = Q_("0.9 kg * m**2")  # Based on slightly different propeller
izg = Q_("0.2 kg * m**2")  # DUMMY

#Propeller Diameter
D_prop = Q_("1.90 m")
