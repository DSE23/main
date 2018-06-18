"""
Name: Airfoil discrete vortex simulation
Department: Aerodynamics
Last updated: 11/06/2018 12:45 by Sam
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

N = 40                         # Number of panels to be used in total
alpha = m.radians(-2.862405226)
Vspeed = 60
Vvector = Vspeed*np.array([m.cos(alpha),0,m.sin(alpha)])
unitx = np.array([1,0,0])
unity = np.array([0,1,0])
unitz = np.array([0,0,1])

def airfoilcamberline(x):
    return 0*x

def transformaileron(xaileron,yaileron,zaileron,da,hingeline):
    points = np.array([xaileron, yaileron, zaileron]).transpose()
    transformationmatrix = np.array([[m.cos(da),0,m.sin(da)],
                                     [0,1,0],
                                     [-m.sin(da),0,m.cos(da)]])
    for i in range(np.size(points,0)):
        points[i] = transformationmatrix.dot(points[i]-hingeline)+hingeline
    return points

def matrixconstructor(normallist,vortexlist,controllist):
    matrix = np.empty((N,N))
    for i in range(N):
        for j in range(N):
            r = controllist[i]-vortexlist[j]
            rnorm = np.linalg.norm(r)
            matrix[i,j] = ((np.cross(-unity,r)/(2*m.pi))/rnorm**2).dot(normallist[i])
    return matrix

def righthandsideconstructor(normallist):
    RHS = np.empty(N)
    for i in range(N):
        RHS[i] = -Vvector.dot(normallist[i])
    return RHS

def computevelocityupper(tangentlist,vortexlist,controlpoint,gammalist):
    V = Vspeed*np.array([m.cos(alpha),0,m.sin(alpha)])
    for i in range(N):
        r = controlpoint - vortexlist[i] + np.array([0,0,0])
        rnorm = np.linalg.norm(r)
        V1 = ((np.cross(-unity, r) / (2 * m.pi)) / rnorm ** 2)
        print(V1)
        V += V1*gammalist[i]
    return V

def computevelocitylower(tangentlist,vortexlist,controlpoint,gammalist):
    V = Vspeed*np.array([m.cos(alpha),0,m.sin(alpha)])
    for i in range(N):
        r = controlpoint - vortexlist[i] - np.array([0,0,0])
        rnorm = np.linalg.norm(r)
        V1 = ((np.cross(-unity, r) / (2 * m.pi)) / rnorm ** 2).dot(tangentlist[i])
        V += V1*gammalist[i]
    return V

def pressurecoefficientuppercalculator(tangentlist,vortexlist,controlpoint,gammalist):
    Vlocal = computevelocityupper(tangentlist,vortexlist,controlpoint,gammalist)
    Vmag = np.linalg.norm(Vlocal)
    return Vlocal

def pressurecoefficientlowercalculator(tangentlist,vortexlist,controlpoint,gammalist):
    Vlocal = computevelocitylower(tangentlist,vortexlist,controlpoint,gammalist)
    Vmag = np.linalg.norm(Vlocal)
    return Vlocal



### Aileron geometry
aileronlength = 0.25
aileroncenter = (2-aileronlength)/2
Naileron = int(aileronlength*N)
hingelocation = 1.0
xhingeline = aileroncenter + (hingelocation-0.5) * aileronlength
da = m.radians(0)

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
panelarray[Nairfoil:N,2,:] = panelarray[Nairfoil:N,1,:] - panelarray[Nairfoil:N,0,:]
panelarray[Nairfoil:N,3,:] = panelarray[Nairfoil:N,2,:]/(np.sqrt((panelarray[Nairfoil:N,2,:]*panelarray[Nairfoil:N,2,:]).sum(axis=1)).reshape(-1,1))
panelarray[Nairfoil:N,4,:] = np.cross(panelarray[Nairfoil:N,3,:],unity)
panelarray[Nairfoil:N,5,:] = 0.25*panelarray[Nairfoil:N,0,:] + 0.75*panelarray[Nairfoil:N,1,:]      #Control points
panelarray[Nairfoil:N,6,:] = 0.75*panelarray[Nairfoil:N,0,:] + 0.25*panelarray[Nairfoil:N,1,:]      #Vortex points


A = matrixconstructor(panelarray[:,4,:],panelarray[:,6,:],panelarray[:,5,:])
# print(A)
RHS = righthandsideconstructor(panelarray[:,4,:])

Gammalist = np.linalg.solve(A,RHS)

Cl = np.sum(Gammalist)/(0.5*Vspeed)
print(Cl)
print(Gammalist)

cpressureupperlist = np.zeros((N,3))
cpressurelowerlist = np.zeros((N,3))


#
for j in range(N):
    cpressureupperlist[j] = pressurecoefficientuppercalculator(panelarray[:,3,:],panelarray[:,6,:],panelarray[j,5,:],Gammalist)
    cpressurelowerlist[j] = pressurecoefficientlowercalculator(panelarray[:,3,:],panelarray[:,6,:],panelarray[j,5,:],Gammalist)



# plt.plot(panelarray[:,6,0],cpressureupperlist,'x',panelarray[:,6,0],cpressurelowerlist,'x')
# plt.show()