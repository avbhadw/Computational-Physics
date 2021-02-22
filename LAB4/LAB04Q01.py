from SolveLinear import GaussElim,PartialPivot
from numpy import*
import matplotlib.pyplot as plt
from numpy.linalg import solve
from numpy.random import rand
from time import time
#define values of N
N = arange(5,201,1)
err_p = []
err_g = []
err_lu = []
timep = []
timeg = []
timelu = []
for n in N:
    A,v = rand(n, n), rand(n)
    #call the functions and time them
    startp = time()
    x_pp= (PartialPivot(A, v))
    endp = time()
    timep.append(endp-startp)
    A,v = rand(n, n), rand(n)
    startg = time()
    x_ge = (GaussElim(A, v))
    endg = time()
    timeg.append(endg - startg)
    A,v = rand(n, n), rand(n)
    startl = time()
    x_lu = (solve(A, v))
    endl = time()
    timelu.append(endl-startl)
    #comparing with true value
    v_sol_p = (dot(A, x_pp))
    v_sol_g = (dot(A, x_ge))
    v_sol_lu = (dot(A, x_lu))
    err_p.append((mean(abs(v-v_sol_p))))
    err_g.append(mean(abs(v-v_sol_g)))
    err_lu.append(mean(abs(v-v_sol_lu)))
plt.plot(N,log(err_g), label='Gaussian')
plt.plot(N,log(err_p), label = 'Partial Pivot')
plt.plot(N,log(err_lu), label = 'LU')
plt.title('Log plot of Relative error vs N ')
plt.xlabel('N')
plt.ylabel('Error')
plt.legend()
plt.show()
plt.savefig('NvslogE')
plt.plot(N,log(timeg), label='Gaussian')
plt.plot(N,log(timep), label = 'Partial Pivot')
plt.plot(N,log(timelu), label = 'LU')
plt.title('Log plot of Time vs N')
plt.xlabel('N')
plt.ylabel('Log(Time(s))')
plt.legend()
plt.show()
plt.savefig('NvslogT')

