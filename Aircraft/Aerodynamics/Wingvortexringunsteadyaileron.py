"""
Name: Airfoil discrete vortex simulation
Department: Aerodynamics
Last updated: 11/06/2018 12:45 by Sam
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()

### Mesh input parameters
Nspanwise = 40
Nchordwise = 20
dt = 0.01

### C.g. definition
positioncg = np.array([0.,0.,0.])

### Input wing geometry
fuselagediameter = 0.
wingspan = 7.54
rootchord = 1.8
tipchord = 0.8
wingzerosweepline = 0.25

### Input aileron geometry
rootpercentage = 0.25
tippercentage = 0.25
ailerondeflection = m.radians(25)
hornwidth = 0.00001
hingelinepercentageroot = 0.8
hingelinepercentagetip = 0.8

### Output geometry
fuselagehalfidameter = fuselagediameter/2
winghalfspan = wingspan/2
effectivewinghalfspan = winghalfspan-fuselagehalfidameter-hornwidth
effectivewinghalfspanincludinghorn = effectivewinghalfspan + hornwidth
winggeometriccenter = fuselagehalfidameter+effectivewinghalfspan/2

### Mesh computations
Nspanwisewing = int((effectivewinghalfspan)/effectivewinghalfspanincludinghorn*Nspanwise)
Nspanwisehorn = Nspanwise-Nspanwisewing
Nchordwisewing = int((1-tippercentage)*Nchordwise)
Nchordwiseaileron = Nchordwise-Nchordwisewing


### Mesh creation
if hornwidth < 0.00001:
    thetaspanwisewing = np.linspace(0, m.pi, Nspanwisewing, endpoint=True)
    thetachordwisewing = np.linspace(0, m.pi, Nchordwisewing, endpoint=False)
    thetachordwiseaileron = np.linspace(0, m.pi, Nchordwiseaileron, endpoint=True)

    thetaspanwisehorn = np.linspace(0, m.pi, Nspanwisehorn, endpoint=True)
    thetachordwisehorn = np.linspace(0, m.pi, Nchordwise, endpoint=True)

else:
    thetaspanwisewing = np.linspace(0,m.pi-0.00001,Nspanwisewing+1,endpoint=True)
    thetachordwisewing = np.linspace(0,m.pi-0.00001,Nchordwisewing+1,endpoint=True)
    thetachordwiseaileron = np.linspace(0,m.pi,Nchordwiseaileron,endpoint=True)

    thetaspanwisehorn = np.linspace(0,m.pi,Nspanwisehorn,endpoint=True)
    thetachordwisehorn = np.linspace(0,m.pi,Nchordwise,endpoint=True)

xwing = np.linspace(0,0.5-0.0000001,Nchordwisewing+1,endpoint=True)  # 0.25 - np.cos(thetachordwisewing)/4
ywing = np.linspace(0,0.5-0.0000001,Nspanwisewing+1,endpoint=True)  # 0.25 - np.cos(thetaspanwisewing)/4
# ywing = winggeometriccenter - effectivewinghalfspan*np.cos(thetaspanwisewing)/2

xaileron = np.linspace(0.5,1,Nchordwiseaileron,endpoint=True)    # 0.75 - np.cos(thetachordwiseaileron)/4
yaileron = np.linspace(0,0.5-0.0000001,Nspanwisewing+1,endpoint=True) # 0.25 - np.cos(thetaspanwisewing)/4

xhorn = 0.5 - np.cos(thetachordwisehorn)/2
yhorn = 0.75 - np.cos(thetaspanwisehorn)/4

xpositions = np.concatenate((xwing,xaileron))
ypositions = np.concatenate((ywing,yhorn))

def computeaileronpercentage(y):
    return rootpercentage - (rootpercentage-tippercentage)*y

def transform(x,y):
    xarray = 0

def hornpercentage():
    return hornwidth/effectivewinghalfspanincludinghorn

def scaletospanwise(y):
    hornfraction = hornpercentage()
    return hornfraction, y*((y<0.5)*(1-hornfraction)/0.5)+(y>=0.499999999999999)*(1-(1-y)*hornfraction/0.5)

def hingeline(y):
    return hingelinepercentageroot - (hingelinepercentageroot-hingelinepercentagetip)*y

def airfoilcamberline(x):
    return 0*x

def transformaileron(array,hingelinearray):
    transformationmatrix = np.array([[m.cos(ailerondeflection),0,-m.sin(ailerondeflection)],
                                     [0,1,0],
                                     [m.sin(ailerondeflection),0,m.cos(ailerondeflection)]])
    #print(np.einsum('lk,ijk->ijl',transformationmatrix,array-hingelinearray)+hingelinearray)
    return np.einsum('lk,ijk->ijl',transformationmatrix,array-hingelinearray)+hingelinearray

def chordlengthunitspan(y):
    return rootchord-(rootchord-tipchord)*y

def addminussign(y):
    return -y

def elongatespan(y):
    return fuselagehalfidameter + y*effectivewinghalfspan

def scaletoaileronpercentage(x,y):
    hornfraction, y = scaletospanwise(y)
    aileronpercentage = computeaileronpercentage(y)
    scaledx = x[:,np.newaxis]*((x[:,np.newaxis]<0.5)*(1-aileronpercentage)/0.5)+(x[:,np.newaxis]>=0.49999999999)*(1-(1-x[:,np.newaxis])*aileronpercentage/0.5)
    z = airfoilcamberline(scaledx)
    array = np.dstack((scaledx,np.tile(y,(Nchordwise+1,1)),z))
    savedarray = array
    aileronpercentagearray = np.dstack((computeaileronpercentage(array[:,:,1]),array[:,:,1],np.zeros(np.shape(array[:,:,2])))) ### Switch x and y around
    hingelinearray = np.dstack((hingeline(array[:,:,1]),array[:,:,1],np.zeros(np.shape(array[:,:,2])))) ### Switch x and y around
    newarray =  transformaileron(array,hingelinearray)
    indicestotake = (savedarray[:,:,0]<(1-aileronpercentagearray[:,:,0]))*(savedarray[:,:,1]<(1-hornfraction))
    newarray[indicestotake] = savedarray[indicestotake]
    chordlengtharray = chordlengthunitspan(newarray[:,:,1])
    chordlenghtmultiplicationarray = np.dstack((chordlengtharray,np.ones(np.shape(chordlengtharray)),chordlengtharray))
    newarray = (newarray-np.array([wingzerosweepline,0,0]))*chordlenghtmultiplicationarray
    newarray[:,:,1] = -elongatespan(newarray[:,:,1])
    newarray[:,:,0] = -newarray[:,:,0]
    return newarray

def createpanelarray(pointsarray):
    panelarray = np.zeros((2*Nspanwise*Nchordwise,15,3))
    panelarray[Nspanwise*Nchordwise:,0,:] = pointsarray[:-1,:-1,:].reshape(-1,3)
    panelarray[Nspanwise*Nchordwise:,1,:] = pointsarray[1:,:-1,:].reshape(-1,3)
    panelarray[Nspanwise*Nchordwise:,2,:] = pointsarray[1:,1:,:].reshape(-1,3)
    panelarray[Nspanwise*Nchordwise:,3,:] = pointsarray[:-1,1:,:].reshape(-1,3)
    panelarray = np.delete(panelarray,range(Nspanwisewing+Nspanwise*Nchordwise,2*Nspanwise*Nchordwise,Nspanwise),axis=0)
    panelarray = np.delete(panelarray,range((Nspanwise-1)*Nchordwisewing+Nspanwise*Nchordwise,(Nspanwise-1)*Nchordwisewing+Nspanwise*Nchordwise+Nspanwise-1),axis=0)
    symmetryarray = np.copy(panelarray[:Nspanwise*Nchordwise-1:-1,:,:])
    panelarray = np.delete(panelarray,range(Nspanwise+Nchordwise-1),axis=0)
    symmetryarray[:,:,1] = addminussign(symmetryarray[:,:,1])
    tempsymmetry1 = np.copy(symmetryarray[:,2:4,:])
    symmetryarray[:,3:1:-1,:] = np.copy(symmetryarray[:,0:2,:])
    symmetryarray[:,1::-1,:] = np.copy(tempsymmetry1)
    panelarray[:Nspanwise*Nchordwise-Nspanwise-Nchordwise+1,:,:] = symmetryarray
    panelarray[:,4,:] = panelarray[:,1,:] - panelarray[:,0,:]
    panelarray[:,5,:] = panelarray[:,2,:] - panelarray[:,3,:]
    panelarray[:,0,:] = panelarray[:,0,:] + 0.25*panelarray[:,4,:]
    panelarray[:,1,:] = panelarray[:,1,:] + 0.25*panelarray[:,4,:]
    panelarray[:,2,:] = panelarray[:,2,:] + 0.25*panelarray[:,5,:]
    panelarray[:,3,:] = panelarray[:,3,:] + 0.25*panelarray[:,5,:]
    panelarray[:,6,:] = panelarray[:,2,:] - panelarray[:,0,:]
    panelarray[:,7,:] = panelarray[:,3,:] - panelarray[:,1,:]
    panelarray[:,8,:] = np.cross(panelarray[:,6,:],panelarray[:,7,:])
    panelarray[:,9,:] = panelarray[:,8,:]/(np.sqrt((panelarray[:,8,:]*panelarray[:,8,:]).sum(axis=1)).reshape(-1,1))
    panelarray[:,10,:] = (panelarray[:,4,:]+panelarray[:,5,:])/2
    panelarray[:,11,:] = (panelarray[:,0,:]+panelarray[:,3,:])/2
    panelarray[:,12,:] = panelarray[:,11,:] + 0.5*panelarray[:,10,:]
    panelarray[:,13,:] = ((panelarray[:,1,:]-panelarray[:,0,:])+(panelarray[:,2,:]-panelarray[:,3,:]))/2                       #Tangential vector in chordwise direction, tau_j
    panelarray[:,14,:] = ((panelarray[:,0,:]-panelarray[:,3,:])+(panelarray[:,1,:]-panelarray[:,2,:]))/2                       #Tangential vector in spanwise direction, tau_i
    panelarray[:,13,:] = panelarray[:,13,:]/np.linalg.norm(panelarray[:,13,:],axis=-1)[:,np.newaxis]
    panelarray[:,14,:] = panelarray[:,14,:]/np.linalg.norm(panelarray[:,14,:],axis=-1)[:,np.newaxis]
    areasofpanels = (np.linalg.norm(np.cross(panelarray[:,1,:]-panelarray[:,0,:],panelarray[:,3,:]-panelarray[:,0,:]),axis=-1)
                          +np.linalg.norm(np.cross(panelarray[:,1,:]-panelarray[:,2,:],panelarray[:,3,:]-panelarray[:,2,:]),axis=-1))/2
    chordsofpanels = ((panelarray[:,1,0]-panelarray[:,0,0])+(panelarray[:,2,0]-panelarray[:,3,0]))/2
    spansofpanels = ((panelarray[:,3,1]-panelarray[:,0,1])+(panelarray[:,2,1]-panelarray[:,1,1]))/2
    panelarray[:Nspanwise-1,1,0] -= 1000000
    panelarray[:Nspanwise-1,2,0] -= 1000000
    panelarray[-Nspanwise+1:,1,0] -= 1000000
    panelarray[-Nspanwise+1:,2,0] -= 1000000
    return panelarray, areasofpanels, chordsofpanels, spansofpanels

pointsarray = scaletoaileronpercentage(xpositions,ypositions)
panelarray, areasofpanels, chordsofpanels, spansofpanels = createpanelarray(pointsarray)

numberofpanelsspanwise = 2*(Nspanwise-1)
numberofpanelsschordwise = Nchordwise-1
listofpanelpoints = panelarray[:,0:4,:]
listofcontrolpoints = panelarray[:,12,:]

a = np.zeros((numberofpanelsschordwise*numberofpanelsspanwise,numberofpanelsschordwise*numberofpanelsspanwise,4,3))

matrixofdistancesofcontrolpoints = listofcontrolpoints[:,np.newaxis,np.newaxis,:] - listofpanelpoints[np.newaxis,:,:,:]
matrixofdistancesofpanelpoints = np.zeros((numberofpanelsschordwise*numberofpanelsspanwise,4,3))
matrixofdistancesofpanelpoints[:,0,:] = panelarray[:,1,:] - panelarray[:,0,:]
matrixofdistancesofpanelpoints[:,1,:] = panelarray[:,2,:] - panelarray[:,1,:]
matrixofdistancesofpanelpoints[:,2,:] = panelarray[:,3,:] - panelarray[:,2,:]
matrixofdistancesofpanelpoints[:,3,:] = panelarray[:,0,:] - panelarray[:,3,:]
# matrixofnormsofdistancesofpanelpoints = np.sqrt((matrixofdistancesofpanelpoints*matrixofdistancesofpanelpoints).sum(axis=2))
matrixofnormsofdistancesofcontrolpoints = np.sqrt((matrixofdistancesofcontrolpoints*matrixofdistancesofcontrolpoints).sum(axis=3))

lefthandsidematrix = 1/(4*m.pi)*np.einsum('ijk,jk->ij',np.cross(matrixofdistancesofcontrolpoints[:,:,0,:],matrixofdistancesofcontrolpoints[:,:,1,:]) / \
                                          (np.linalg.norm(np.cross(matrixofdistancesofcontrolpoints[:, :, 0, :],
                                                                  matrixofdistancesofcontrolpoints[:, :, 1, :]),
                                                         axis=-1)**2)[:,:,np.newaxis] * \
                     np.einsum('jk,ijk->ij',matrixofdistancesofpanelpoints[:,0,:],matrixofdistancesofcontrolpoints[:,:,0,:]/
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,0,np.newaxis]-
                                                                           matrixofdistancesofcontrolpoints[:,:,1,:] /
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,1,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])
lefthandsidematrix += 1/(4*m.pi)*np.einsum('ijk,jk->ij',np.cross(matrixofdistancesofcontrolpoints[:,:,1,:],matrixofdistancesofcontrolpoints[:,:,2,:]) / \
                                          (np.linalg.norm(np.cross(matrixofdistancesofcontrolpoints[:, :, 1, :],
                                                                  matrixofdistancesofcontrolpoints[:, :, 2, :]),
                                                         axis=-1)**2)[:,:,np.newaxis] * \
                     np.einsum('jk,ijk->ij',matrixofdistancesofpanelpoints[:,1,:],matrixofdistancesofcontrolpoints[:,:,1,:]/
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,1,np.newaxis]-
                                                                           matrixofdistancesofcontrolpoints[:,:,2,:] /
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,2,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])

lefthandsidematrix += 1/(4*m.pi)*np.einsum('ijk,jk->ij',np.cross(matrixofdistancesofcontrolpoints[:,:,2,:],matrixofdistancesofcontrolpoints[:,:,3,:]) / \
                                          (np.linalg.norm(np.cross(matrixofdistancesofcontrolpoints[:, :, 2, :],
                                                                  matrixofdistancesofcontrolpoints[:, :, 3, :]),
                                                         axis=-1)**2)[:,:,np.newaxis] * \
                     np.einsum('jk,ijk->ij',matrixofdistancesofpanelpoints[:,2,:],matrixofdistancesofcontrolpoints[:,:,2,:]/
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,2,np.newaxis]-
                                                                           matrixofdistancesofcontrolpoints[:,:,3,:] /
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,3,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])

lefthandsidematrix += 1/(4*m.pi)*np.einsum('ijk,jk->ij',np.cross(matrixofdistancesofcontrolpoints[:,:,3,:],matrixofdistancesofcontrolpoints[:,:,0,:]) / \
                                          (np.linalg.norm(np.cross(matrixofdistancesofcontrolpoints[:, :, 3, :],
                                                                  matrixofdistancesofcontrolpoints[:, :, 0, :]),
                                                         axis=-1)**2)[:,:,np.newaxis] * \
                     np.einsum('jk,ijk->ij',matrixofdistancesofpanelpoints[:,3,:],matrixofdistancesofcontrolpoints[:,:,3,:]/
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,3,np.newaxis]-
                                                                           matrixofdistancesofcontrolpoints[:,:,0,:] /
                                                                           matrixofnormsofdistancesofcontrolpoints[:,:,0,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])

wakestrengths = np.zeros((numberofpanelsspanwise*20))
wakepanels = np.random.rand(numberofpanelsspanwise*20,4,3)
cgvelocity = np.array([60,0,4])
Omega = np.array([0,0,0])
velocityduetoangularrates = np.cross(Omega,panelarray[:,12,:])
localvelocity = cgvelocity+velocityduetoangularrates
righthandsidevector = np.einsum('ij,ij->i',localvelocity,panelarray[:,9,:])

solution = np.linalg.solve(lefthandsidematrix,righthandsidevector)
t = 0
while t < 0.01:

    wakestrengths[:-numberofpanelsspanwise] = wakestrengths[numberofpanelsspanwise:]
    wakestrengths[-numberofpanelsspanwise:] = np.concatenate((solution[:numberofpanelsspanwise//2],solution[-numberofpanelsspanwise//2:]))
    wakepanels[:-numberofpanelsspanwise,:,:] = wakepanels[numberofpanelsspanwise:,:,:]

    trailingedgepanels = np.concatenate((panelarray[:numberofpanelsspanwise//2,:,:],panelarray[-numberofpanelsspanwise//2:,:,:]),axis=0)

    wakepanels[-numberofpanelsspanwise:,0,:] = trailingedgepanels[:,1,:]
    # wakepanels[:,1,:] = panelarray[:numberofpanelsschordwise,1,:]
    # wakepanels[:,2,:] = panelarray[:numberofpanelsschordwise,2,:]
    wakepanels[-numberofpanelsspanwise:,3,:] = trailingedgepanels[:,2,:]

    panelarray[:,0:4,:] += cgvelocity[np.newaxis,np.newaxis,:]*dt + np.cross(Omega,panelarray[:,0:4,:]-positioncg[np.newaxis,np.newaxis,:])
    panelarray[:,12,:] += cgvelocity[np.newaxis,:]*dt + np.cross(Omega,panelarray[:,12,:]-positioncg[np.newaxis,:])
    trailingedgepanels = np.concatenate((panelarray[:numberofpanelsspanwise//2,:,:],panelarray[-numberofpanelsspanwise//2:,:,:]),axis=0)
    wakepanels[-numberofpanelsspanwise:,1,:] = trailingedgepanels[:,1,:]
    wakepanels[-numberofpanelsspanwise:,2,:] = trailingedgepanels[:,2,:]

    velocityduetoangularrates = np.cross(Omega,(panelarray[:,12,:]-positioncg))
    localvelocity = cgvelocity+velocityduetoangularrates
    righthandsidevector = np.einsum('ij,ij->i',localvelocity,panelarray[:,9,:])

    matrixofdistancesofcontrolpointstowakepoints = listofcontrolpoints[:,np.newaxis,np.newaxis,:] - wakepanels[np.newaxis,:,:,:]
    matrixofdistancesofwakepoints = np.zeros((20*numberofpanelsspanwise,4,3))
    matrixofdistancesofwakepoints[:,0,:] = wakepanels[:,1,:] - wakepanels[:,0,:]
    matrixofdistancesofwakepoints[:,1,:] = wakepanels[:,2,:] - wakepanels[:,1,:]
    matrixofdistancesofwakepoints[:,2,:] = wakepanels[:,3,:] - wakepanels[:,2,:]
    matrixofdistancesofwakepoints[:,3,:] = wakepanels[:,0,:] - wakepanels[:,3,:]
    matrixofnormsofdistancesofcontrolpointstowakepoints = np.sqrt((matrixofdistancesofcontrolpointstowakepoints*matrixofdistancesofcontrolpointstowakepoints).sum(axis=3))

    matrixofinducedvelocityofwake = 1/(4*m.pi)*np.einsum('ijk,ik->ij',np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,0,:],matrixofdistancesofcontrolpointstowakepoints[:,:,1,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 0, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 1, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                            np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,0,:],matrixofdistancesofcontrolpointstowakepoints[:,:,0,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,0,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,1,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,1,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])
    matrixofinducedvelocityofwake += 1/(4*m.pi)*np.einsum('ijk,ik->ij',np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,1,:],matrixofdistancesofcontrolpointstowakepoints[:,:,2,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 1, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 2, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                         np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,1,:],matrixofdistancesofcontrolpointstowakepoints[:,:,1,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,1,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,2,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,2,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])

    matrixofinducedvelocityofwake += 1/(4*m.pi)*np.einsum('ijk,ik->ij',np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,2,:],matrixofdistancesofcontrolpointstowakepoints[:,:,3,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 2, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 3, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                         np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,2,:],matrixofdistancesofcontrolpointstowakepoints[:,:,2,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,2,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,3,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,3,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])

    matrixofinducedvelocityofwake += 1/(4*m.pi)*np.einsum('ijk,ik->ij',np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,3,:],matrixofdistancesofcontrolpointstowakepoints[:,:,0,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 3, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 0, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                         np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,3,:],matrixofdistancesofcontrolpointstowakepoints[:,:,3,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,3,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,0,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,0,np.newaxis])[:,:,np.newaxis],panelarray[:,9,:])

    righthandsidevector -= np.einsum('ij,j->i',matrixofinducedvelocityofwake,wakestrengths)
    solutionnew = np.linalg.solve(lefthandsidematrix,righthandsidevector)

    print("Error after", t, "is",np.linalg.norm(solutionnew-solution))

    matrixofinducedvelocitywakereadyfortangential = 1/(4*m.pi)*np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,0,:],matrixofdistancesofcontrolpointstowakepoints[:,:,1,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 0, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 1, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                            np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,0,:],matrixofdistancesofcontrolpointstowakepoints[:,:,0,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,0,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,1,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,1,np.newaxis])[:,:,np.newaxis]

    matrixofinducedvelocitywakereadyfortangential += 1/(4*m.pi)*np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,1,:],matrixofdistancesofcontrolpointstowakepoints[:,:,2,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 1, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 2, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                         np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,1,:],matrixofdistancesofcontrolpointstowakepoints[:,:,1,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,1,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,2,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,2,np.newaxis])[:,:,np.newaxis]

    matrixofinducedvelocitywakereadyfortangential += 1/(4*m.pi)*np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,2,:],matrixofdistancesofcontrolpointstowakepoints[:,:,3,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 2, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 3, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                         np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,2,:],matrixofdistancesofcontrolpointstowakepoints[:,:,2,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,2,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,3,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,3,np.newaxis])[:,:,np.newaxis]

    matrixofinducedvelocitywakereadyfortangential += 1/(4*m.pi)*np.cross(matrixofdistancesofcontrolpointstowakepoints[:,:,3,:],matrixofdistancesofcontrolpointstowakepoints[:,:,0,:]) / \
                                              (np.linalg.norm(np.cross(matrixofdistancesofcontrolpointstowakepoints[:, :, 3, :],
                                                                       matrixofdistancesofcontrolpointstowakepoints[:, :, 0, :]),
                                                             axis=-1)**2)[:,:,np.newaxis] * \
                         np.einsum('jk,ijk->ij',matrixofdistancesofwakepoints[:,3,:],matrixofdistancesofcontrolpointstowakepoints[:,:,3,:]/
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,3,np.newaxis]-
                                   matrixofdistancesofcontrolpointstowakepoints[:,:,0,:] /
                                   matrixofnormsofdistancesofcontrolpointstowakepoints[:,:,0,np.newaxis])[:,:,np.newaxis]
    matrixofinducedvelocitywakereadyfortangential = np.einsum('ijk,j->ik',matrixofinducedvelocitywakereadyfortangential,wakestrengths)+localvelocity

    gammadifferenceinchordwise = np.zeros(numberofpanelsspanwise*numberofpanelsschordwise)
    gammadifferenceinchordwise[0:numberofpanelsspanwise//2*(numberofpanelsschordwise-1)] = solution[0:numberofpanelsspanwise//2*(numberofpanelsschordwise-1)] \
                                                                                        -solution[numberofpanelsspanwise//2:numberofpanelsspanwise//2*numberofpanelsschordwise]
    gammadifferenceinchordwise[numberofpanelsspanwise//2*(numberofpanelsschordwise-1):numberofpanelsspanwise//2*numberofpanelsschordwise] \
        = solution[numberofpanelsspanwise//2*(numberofpanelsschordwise-1):numberofpanelsspanwise//2*numberofpanelsschordwise]
    gammadifferenceinchordwise[numberofpanelsspanwise//2*(numberofpanelsschordwise):numberofpanelsspanwise//2*(numberofpanelsschordwise+1)] \
        = solution[numberofpanelsspanwise//2*(numberofpanelsschordwise):numberofpanelsspanwise//2*(numberofpanelsschordwise+1)]
    gammadifferenceinchordwise[numberofpanelsspanwise//2*(numberofpanelsschordwise+1):] = solution[numberofpanelsspanwise//2*(numberofpanelsschordwise+1):]\
                                                                                          -solution[numberofpanelsspanwise//2*numberofpanelsschordwise:-numberofpanelsspanwise//2]

    gammadifferenceinspanwise = np.zeros(numberofpanelsspanwise*numberofpanelsschordwise)
    gammadifferenceinspanwise[1:numberofpanelsspanwise//2*numberofpanelsschordwise] = solution[:numberofpanelsspanwise//2*numberofpanelsschordwise-1]\
                                                                                      -solution[1:numberofpanelsspanwise//2*numberofpanelsschordwise]
    gammadifferenceinspanwise[0:numberofpanelsspanwise//2*numberofpanelsschordwise:numberofpanelsspanwise//2] \
        = solution[0:numberofpanelsspanwise//2*numberofpanelsschordwise:numberofpanelsspanwise//2]\
          -solution[1:numberofpanelsspanwise//2*numberofpanelsschordwise:numberofpanelsspanwise//2]
    gammadifferenceinspanwise[numberofpanelsspanwise//2*numberofpanelsschordwise:-1] = solution[numberofpanelsspanwise//2*numberofpanelsschordwise:-1]\
                                                                                       -solution[numberofpanelsspanwise//2*numberofpanelsschordwise+1:]
    gammadifferenceinspanwise[numberofpanelsspanwise//2*(numberofpanelsschordwise+1)-1::numberofpanelsspanwise//2] = \
        solution[numberofpanelsspanwise//2*(numberofpanelsschordwise+1)-2::numberofpanelsspanwise//2]-\
        solution[numberofpanelsspanwise//2*(numberofpanelsschordwise+1)-1::numberofpanelsspanwise//2]

    pressure = 1.225*(np.einsum('ik,ik->i',matrixofinducedvelocitywakereadyfortangential,panelarray[:,13,:])*gammadifferenceinchordwise/chordsofpanels+
                      np.einsum('ik,ik->i', matrixofinducedvelocitywakereadyfortangential,
                                panelarray[:, 14, :])*gammadifferenceinspanwise/spansofpanels+(solutionnew-solution)*dt)
    lift = 1.225*60*gammadifferenceinchordwise*spansofpanels
    print(np.sum(lift))
    F = (pressure*areasofpanels)[:,np.newaxis]*panelarray[:,9,:]
    F = np.sum(F,axis=0)
    print(F)
    solution = solutionnew
    t+=dt




# lift = (1.225*60*solution*(panelarray[:,2,1]-panelarray[:,0,1])/\
#       (np.linalg.norm(np.cross(panelarray[:,1,:]-panelarray[:,0,:],panelarray[:,3,:]-panelarray[:,0,:]),axis=-1)
#        +np.linalg.norm(np.cross(panelarray[:,1,:]-panelarray[:,2,:],panelarray[:,3,:]-panelarray[:,2,:]),axis=-1))/2)

end = time.time()
print(end-start)