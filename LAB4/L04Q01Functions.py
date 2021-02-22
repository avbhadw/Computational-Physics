import numpy as np
import matplotlib.pyplot  as plt
#solving linear systems
import numpy as np
from time import time
# A1 = np.array([[ 2,  1,  4,  1 ],
#            [ 3,  4, -1, -1 ],
#            [ 1, -4,  1,  5 ],
#            [ 2, -2,  1,  3 ]],float)
# v1 = np.array([ -4, 3, 9, 7 ],float)

# Gaussian elimination
def GaussianElim(A,v):
    N = len(v)
    for m in range(N):

    # Divide by the diagonal element
        div = A[m,m]
        A[m,:] /= div
        v[m] /= div

    # Now subtract from the lower rows
        for i in range(m+1,N):
            mult = A[i,m]
            A[i,:] -= mult*A[m,:]
            v[i] -= mult*v[m]

# Backsubstitution
    x = np.empty(N,float)
    for m in range(N-1,-1,-1):
        x[m] = v[m]
        for i in range(m+1,N):
            x[m] -= A[m,i]*x[i]
    return x

def PartialPivot(A,v):
    N= len(v)
    for m in range (len(v)):
        for j in range(m,len(v)):
            if A[j,m]>A[m,m]:
                #swap matrices here
                A[m,:],A[:,j] = np.copy(A[m,:]),np.copy(A[j,:])
                v[m,:],v[:,j] = np.copy(v[j,:]),np.copy[m,:]
            # Divide by the diagonal element
            div = A[m, m]
            A[m, :] /= div
            v[m] /= div
            # Now subtract from the lower rows
            for i in range(m + 1, N):
                mult = A[i, m]
                A[i, :] -= mult * A[m, :]
                v[i] -= mult * v[m]
    # Backsubstitution
    x = np.empty(N, complex)
    for m in range(N - 1, -1, -1):
        x[m] = v[m]
        for i in range(m + 1, N):
            x[m] -= A[m, i] * x[i]

    return x

# test = PartialPivot(A1,v1)
# vsol = np.dot(A1,test)
# err = (np.mean(abs(v1-vsol)))
# test2 = GaussianElim(A1,v1)
# # print(err)
# print(test)
# test2 = GaussianElim(A1,v1)
# print(test2)
# # def LU(A,v):

