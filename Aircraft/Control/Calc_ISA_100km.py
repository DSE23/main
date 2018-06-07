from math import *

def isacal(altitude):
    T0 = 288.15    #K
    p0 = 101325.0  #Pa
    rho0 = 1.225   #kg/m3
    g0 = 9.80665   #m/s2
    R = 287.0      #J/kgK
    ft = .3048     #

    T11 = T0  - .0065 * 11000
    T20 = T11
    T32 = T20 + .001  * 12000
    T47 = T32 + .0028 * 15000
    T51 = T47
    T71 = T51 - .0028 * 20000
    T84 = T71 - .002  * 13852

    p11 = (T11/T0)**(-g0/(-.0065*R))*p0
    p20 = exp(-g0/(R*T11)*(20000-11000))*p11
    p32 = (T32/T11)**(-g0/(0.0010*R))*p20
    p47 = (T47/T32)**(-g0/(0.0028*R))*p32
    p51 = exp(-g0/(R*T11)*(51000-47000))*p47
    p71 = (T71/T51)**(-g0/(-.0028*R))*p51
    p84 = (T84/T71)**(-g0/(-.0020*R))*p71

    rho11 = rho0*(T11/T0)**(-g0/(-.0065*R)-1)
    rho20 = exp(-g0/(R*T11)*(20000-11000))*rho11
    rho32 = (T32/T20)**(-g0/(0.0010*R)-1)*rho20
    rho47 = (T47/T32)**(-g0/(0.0028*R)-1)*rho32
    rho51 = exp(-g0/(R*T47)*(51000-47000))*rho47
    rho71 = (T71/T51)**(-g0/(-.0028*R)-1)*rho51
    rho84 = (T84/T71)**(-g0/(-.0020*R)-1)*rho71

    if -610. <= altitude < 11000:
        T1   = T0 - .0065* altitude
        p1   = (T1/T0)**(-g0/(-.0065*R))*p0
        rho1 = ((T1/T0)**(((-g0)/(-.0065*R))-1))*rho0
            
    elif 11000. <= altitude <= 20000.:
        T1   = T11
        p1   = exp(-g0/(R*T11)*(altitude-11000))*p11 
        rho1 = exp(-g0/(R*T11)*(altitude-11000))*rho11
        
    elif 20000. <= altitude < 32000.:
        T1   = T20 + .001 * (altitude-20000)
        p1   = (T1/T20)**(-g0/(.0010*R))*p20
        rho1 = rho20*(T1/T20)**(-g0/(.001*R)-1)
    elif 32000. <= altitude < 47000.:
        T1   = T32 + .0028 * (altitude - 32000)
        p1   = (T1/T32)**(-g0/(.0028*R))*p32
        rho1 = ((T1/T32)**(((-g0)/(.0028*R))-1))*rho32
        
    elif 47000. <= altitude < 51000.:
        T1   = T47
        p1   = exp(-g0/(R*T47)*(altitude-47000))*p47
        rho1 = exp(-g0/(R*T47)*(altitude-47000))*rho47
        
    elif 51000. <= altitude < 71000.:
        T1   = T51 - .0028 * (altitude - 51000)
        p1   = (T1/T51)**(-g0/(-.0028*R))*p51
        rho1 = ((T1/T51)**(((-g0)/(-.0028*R))-1))*rho51
        
    elif 71000. <= altitude < 84852.:
        T1   = T71 - .002 * (altitude - 71000)
        p1   = (T1/T71)**(-g0/(-.0020*R))*p71
        rho1 = ((T1/T71)**(((-g0)/(-.0020*R))-1))*rho71
        
    elif 84852. <= altitude < 100000.:
        T1   = T84
        p1   = exp(-g0/(R*T84)*(altitude-84852))*p84
        rho1 = exp(-g0/(R*T84)*(altitude-84852))*rho84
        
    else:
        T1   = -10
        p1   = -10
        rho1 = -10

    return(T1,p1,rho1)
