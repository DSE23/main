import numpy as np
import time

A = np.array([[0.5,1.5],[2.5,3.5]])
B = np.array([[0,1],[2,3]])
C = A[:,np.newaxis, :] - B[np.newaxis, :, : ]

print(C)
# v = np.arange(0,N)

# M = np.arange(1,26).reshape(5,5)
# A = np.arange(4).reshape((2,2))
# B = np.arange(10,14).reshape((2,2))
#
# s = np.einsum("a->",    v)
# T = np.einsum("ij->ji", M)
# C = np.einsum("mn,np->mp",A,B)