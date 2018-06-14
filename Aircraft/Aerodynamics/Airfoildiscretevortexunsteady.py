"""
Name: Airfoil discrete vortex simulation
Department: Aerodynamics
Last updated: 11/06/2018 12:45 by Sam
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

### Geometry
aileronlength = 0.25
hingelocation = 1.0

### Performance parameters
alpha = m.radians(1)
Vspeed = 60
da = m.radians(30)

### Simulation parameters
N = 40                         # Number of panels to be used in total
Nwakes = 10
t = 0
dt = 0.01
tmax = 20

### Basics
Vvector = Vspeed*np.array([m.cos(alpha),0,m.sin(alpha)])
unitx = np.array([1,0,0])
unity = np.array([0,1,0])
unitz = np.array([0,0,1])

### Airfoil camberline
def airfoilcamberline(x):
    return 0*x

### Transform the coordinates of the aileron based on the location of the hingeline and deflection
def transformaileron(xaileron,yaileron,zaileron,da,hingeline):
    points = np.array([xaileron, yaileron, zaileron]).transpose()
    transformationmatrix = np.array([[m.cos(da),0,m.sin(da)],
                                     [0,1,0],
                                     [-m.sin(da),0,m.cos(da)]])
    for i in range(np.size(points,0)):
        points[i] = transformationmatrix.dot(points[i]-hingeline)+hingeline
    return points

### Function to transform an arbitrary point in the local CS to the global CS
def transformlocaltoglobal(theta,point,origin):
    transformationmatrix = np.array([[m.cos(theta), 0, m.sin(theta)],
                                     [0,1,0],
                                     [-m.sin(theta), 0, m.cos(theta)]])
    return transformationmatrix.dot(point) + origin

### Function to transform an arbitrary point in the global CS to a local CS
def transformglobaltolocal(theta,point,origin):
    transformationmatrix = np.array([[m.cos(theta), 0, -m.sin(theta)],
                                     [0,1,0],
                                     [m.sin(theta), 0, m.cos(theta)]])
    return transformationmatrix.dot(point-origin)

### Compute the matrix coefficients
def matrixcoefficientcalculator(vortexpoint,controlpoint,normal):
    r = controlpoint - vortexpoint
    rnorm = np.linalg.norm(r)
    return ((np.cross(unity, r) / (2 * m.pi)) / rnorm ** 2).dot(normal)

### Construct the matrix
def matrixconstructor(vortexlist,controllist,normallist,lastwakepoint):
    vortexlist = np.append(vortexlist,np.array([[0,0,0]]),axis = 0)
    controllist = np.append(controllist,[lastwakepoint],axis=0)
    normallist = np.append(normallist,np.array([[0,0,0]]),axis =0)
    matrixofdistances = controllist[np.newaxis, :, :] - vortexlist[:, np.newaxis, :]
    normaliseddistances = matrixofdistances / (matrixofdistances * matrixofdistances).sum(axis=2)[:, :, np.newaxis]
    matrixaftercrossproduct = np.cross(unity,normaliseddistances)
    matrixafternormal = 1/(2*m.pi)*np.einsum("ijk,ik->ij",matrixaftercrossproduct,normallist)
    matrixafternormal[-1,:] = 1
    return matrixafternormal

    # matrix = np.ones((N+1,N+1))
    # for i in range(N):
    #     for j in range(N):
    #         matrix[i,j] = matrixcoefficientcalculator(vortexlist[j],controllist[i],normallist[i])
    #     matrix[i,-1] = matrixcoefficientcalculator(lastwakepoint,controllist[i],normallist[i])

### Construct the righthandside
def righthandsideconstructor(normallist,listofwakes,previousgamma):
    RHS = np.empty(N+1)
    for i in range(N):
        RHS = 1
    RHS[N+1] = previousgamma

### Velocity of airfoil origin, as measured in inertial reference frame
def dotX0(t):
    return -Vspeed - 5*m.sin(t)

def dotZ0(t):
    return 3*m.sin(t)

def dottheta0(t):
    return 0.1*m.cos(t)

### Aileron geometry
aileroncenter = (2-aileronlength)/2
Naileron = int(aileronlength*N)
xhingeline = aileroncenter + (hingelocation-0.5) * aileronlength

### Create panel coordinates for the aileron
thetaaileron = np.linspace(0, m.pi,Naileron+1,endpoint=False)
xaileron = aileroncenter - aileronlength * np.cos(thetaaileron)/2
yaileron = np.zeros(Naileron+1)
zaileron = airfoilcamberline(xaileron)
hingeline = np.array([xhingeline,0,airfoilcamberline(xhingeline)])
points = transformaileron(xaileron,yaileron,zaileron,da,hingeline)
xaileron = points[:,0]
yaileron = points[:,1]
zaileron = points[:,2]

### Airfoil geometry
airfoillength = 1 - aileronlength
airfoilcenter = airfoillength/2
Nairfoil = int(airfoillength*N)

### Create panel coordinates for the airfoil
thetaairfoil = np.linspace(0, m.pi,Nairfoil+1,endpoint=False)
xairfoil = airfoilcenter - airfoillength * np.cos(thetaairfoil)/2
yairfoil = np.zeros(Nairfoil+1)
zairfoil = airfoilcamberline(xairfoil)


panelarray = np.empty((N,7,3))
panelarray[0:Nairfoil,0:2,:] = np.dstack((np.column_stack((xairfoil[0:Nairfoil],xairfoil[1:Nairfoil+1])),
                               np.column_stack((yairfoil[0:Nairfoil],yairfoil[1:Nairfoil+1])),
                               np.column_stack((zairfoil[0:Nairfoil],zairfoil[1:Nairfoil+1]))))
panelarray[0:Nairfoil,2,:] = panelarray[0:Nairfoil,0,:] - panelarray[0:Nairfoil,1,:]
panelarray[0:Nairfoil,3,:] = panelarray[0:Nairfoil,2,:]/(np.sqrt((panelarray[0:Nairfoil,2,:]*panelarray[0:Nairfoil,2,:]).sum(axis=1)).reshape(-1,1))
panelarray[0:Nairfoil,4,:] = np.cross(panelarray[0:Nairfoil,3,:],unity)
panelarray[0:Nairfoil,5,:] = 0.25*panelarray[0:Nairfoil,0,:] + 0.75*panelarray[0:Nairfoil,1,:]
panelarray[0:Nairfoil,6,:] = 0.75*panelarray[0:Nairfoil,0,:] + 0.25*panelarray[0:Nairfoil,1,:]

panelarray[Nairfoil:N,0:2,:] = np.dstack((np.column_stack((xaileron[0:Naileron],xaileron[1:Naileron+1])),
                               np.column_stack((yaileron[0:Naileron],yaileron[1:Naileron+1])),
                               np.column_stack((zaileron[0:Naileron],zaileron[1:Naileron+1]))))
panelarray[Nairfoil:N,2,:] = panelarray[Nairfoil:N,0,:] - panelarray[Nairfoil:N,1,:]
panelarray[Nairfoil:N,3,:] = panelarray[Nairfoil:N,2,:]/(np.sqrt((panelarray[Nairfoil:N,2,:]*panelarray[Nairfoil:N,2,:]).sum(axis=1)).reshape(-1,1))
panelarray[Nairfoil:N,4,:] = np.cross(panelarray[Nairfoil:N,3,:],unity)
panelarray[Nairfoil:N,5,:] = 0.25*panelarray[Nairfoil:N,0,:] + 0.75*panelarray[Nairfoil:N,1,:]
panelarray[Nairfoil:N,6,:] = 0.75*panelarray[Nairfoil:N,0,:] + 0.25*panelarray[Nairfoil:N,1,:]

listofwakes = np.zeros((Nwakes,4))
localcoordinatesofwake = np.zeros((Nwakes,3))
airfoilorigin = np.array([0.,0.,0.])
theta = 0
while t < tmax:
    ### Find coordinates of wake in local coordinate system
    for i in range(Nwakes):
        localcoordinatesofwake[i,:]=transformglobaltolocal(theta,listofwakes[i,1:],airfoilorigin)

    ### Construct matrixcoefficients
    matrix = matrixconstructor(panelarray[:,5,:],panelarray[:,6,:],panelarray[:,4,:],localcoordinatesofwake[-1])
    print(t)

    ### Compute RHS vector

    ### Solve system

    ### Compute pressure distribution

    ### Compute velocity at wake points

    ### Update position of wake in airfoil coordinate system


    ### Create new wake point
    localcoordinatesofwake[:-1] = localcoordinatesofwake[1:]
    Omega = np.array([0,dottheta0(t),0])
    VTE = np.array([dotX0(t),0,dotZ0(t)])+Omega.dot(panelarray[-1,1,:])
    localcoordinatesofwake[-1] = panelarray[-1,1,:] + 0.25*panelarray[-1,3,:]*np.linalg.norm(VTE)*dt

    ### Convert positions of wake to global coordinate system
    for i in range(Nwakes):
        listofwakes[i,1:] = transformlocaltoglobal(theta,localcoordinatesofwake[i,:],airfoilorigin)

    ### Update position of airfoil in global coordinate system
    Vairfoilorigin = np.array([dotX0(t),0,dotZ0(t)])
    airfoilorigin += Vairfoilorigin*dt
    t+=dt