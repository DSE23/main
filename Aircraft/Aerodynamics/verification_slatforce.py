# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 10:21:46 2018

@author: Emma
"""
import sys
sys.path.append('../') # This makes sure the parent directory gets added to the system path
from Misc import ureg, Q_
import numpy as np
import math as m
from Geometry import Geometry


#DATCOM graph interpolation
maxLeffect_datapoints = np.array([(0,0),(.05,.84),(.1,1.2),(.15,1.44),(.2,1.6),\
                                  (.25,1.75),(.3,1.825)])
maxL_x = maxLeffect_datapoints[:,0]
maxL_y = maxLeffect_datapoints[:,1]
eta_max = .99
etadelta_datapoints = np.array([(15,1),(17.5,.95),(20,.9),(22.5,.82),(25,7.6),(30,.59)])
ed_x = etadelta_datapoints[:,0]
ed_y = etadelta_datapoints[:,1]

#optimization loop
defl = 30
maxlifteffect = np.interp(.15,maxL_x,maxL_y)
da_corr_fact = np.interp(defl,ed_x,ed_y)
delta_clmax = maxlifteffect * eta_max * .82 * m.radians(defl)

etadelta_datapoints = np.array([(15,1),(17.5,.95),(20,.9),(22.5,.82),(25,.76),(30,.59)])
ed_x = etadelta_datapoints[:,0]
ed_y = etadelta_datapoints[:,1]
            
            


