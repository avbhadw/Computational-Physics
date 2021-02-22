from numpy import linspace
import matplotlib.pyplot  as plt
from math import exp
#define function and derivative of function
def f(x,c):
    return 1 - exp(-c * x)
def dervf(x,c):
    return c * exp(-c * x)
def relaxation(c,accuracy):
    x1 = 1.0
    error = 1.0
    iterations = 1
    while error>accuracy:
        x1, x2 = f(x1,c), x1
        error = abs((x1 - x2) / (1 - 1 / dervf(x2,c)))
        iterations += 1
    return x1,iterations
def overrelaxation(c,accuracy,w):
    x1 = 1.0
    #w = 0.6 #this is the value that gives half as many iterations
    error = 1
    iterations = 1
    while error > accuracy:
        x1, x2 = (1+w)*f(x1, c)-w*x1, x1
        error = abs((x1 - x2) / (1 - 1 /(1+w)*dervf(x2, c)-w))
        iterations += 1
    return x1,iterations
#part a
#part 6.10a
relax = relaxation(2.0,1e-6)
print("x converges to "+str(relax[0]))
print("the number of iterations it takes is:" + str(relax[1]))
#part 6.10b
c = 3.0
#stores values of x
y = []
#array so that c ranges from 0 to 3
constants = linspace(0.001,c,100)
for i in constants:
    test=relaxation(i,1e-6)
    y.append(test[0])
# Plotting
plt.plot(constants,y)
plt.xlabel("c")
plt.ylabel("x")
plt.title('x as a function of c, for $x = 1-e^{-cx}$')
plt.savefig("3a")
#part b

#6.11 a,b

#experimenting with different values of omega
omega = linspace(0,1,10)
result = [15] #no of iterations from relaxation method
lowest = 0
for v in range(len(omega)):
    ovrelax = overrelaxation(2,1e-6,omega[v])
    #finding the lowest no of iterations
    if result[v] > ovrelax[1]:
        lowest = omega[v]
    result.append(ovrelax[1])

overelax = overrelaxation(2,1e-6,lowest)
print("x converges to "+str(overelax[0]))
print('The number of iterations for x to converge is '+ str(overelax[1]))
