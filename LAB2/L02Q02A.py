import numpy as np
import math
from scipy import special
from time import time
def f(t):
    "Returns the dawson function"
    return math.exp(-(4**2))*(math.exp(t ** 2))
def trapezoidal(N,a,b,T):
    """This function calculates the trapezoidal method
    --------------------
    :parameter
    N: the number of slices, float
    a : Start value, float
    b: End value, float
    T: which values to return in the function, boolean
    --------------------
    :returns
    T = 0: the solution to the integral and the comparison with true value (part 2a-i)
    T = 1: the solution to the integral and the time it takes to complete the integral (part 2a-ii)
    T = 2: the solution integral itself (part 2a-iii)"""

    #calculate the width between two points
    h=(b-a)/N
    #store values for f(a) and f(b)
    s = 0.5 * f(a) + 0.5 * f(b)

    startT= time()

    for k in range(1, N):
        s += f(a + k * h)
    sol = h*s

    endT = time()

    diffT = endT - startT
    #difference between true value and numerical approximation
    compareT= abs(sol - special.dawsn(b))

    if (T == 0):
        return sol, compareT
    if (T == 1):
        return sol, diffT
    if (T ==2):
        return sol
def simpson(N,a,b, T):
    """This function calculates the trapezoidal method
       --------------------
       :parameter
       N: the number of slices, float
       a : Start value, float
       b: End value, float
       T: which values to return in the function, boolean
       --------------------
       :returns
       T = 0: the solution to the integral and the comparison with true value (part 2a-i)
       T = 1: the solution to the integral and the time it takes to complete the integral (part 2a-ii)
       T = 2: the solution integral itself (part 2a-iii)"""

    #width of each slice
    h = (b - a) / N

    even =0
    odd = 0

    start2 = time()
    #for odd values
    for k in range(1,N,2):
        odd += 4 * f(a+k*h)
    #for even values
    for k in range (2,N,2):
        even += 2 * f(a+k*h)
    integral = (1/3*h*(f(a)+f(b)+odd+even))
    end2 = time()
    #stores the time
    diff2 = end2 - start2
    #stores the difference between the true value and numerical approximation
    compareS= abs(integral - special.dawsn(b))
    if (T==0) :
        return integral,compareS
    if (T==1):
        return integral, diff2
    if (T==2):
        return integral


#Question 2 part a-i
simprule =simpson(8,0,4,0)
print("Simpson's method with 8 slices  "+ str(simprule[0]))
print("Difference between true value and Simpson's Method " + str(simprule[1]))
trape = trapezoidal(8,0,4,0)
print("Trapezoidal Method with 8 slices"+ str(trape[0]))
print("Difference between true value and Trapezoidal Method " + str(trape[1]))

#2 part a-ii
def ordererror(oerror):
    """The function returns the number of slices needed to have an error within a given order
     --------------------
     :parameter
     oerror: The order of error, float
    --------------------
     :returns
     The number of slices in log2
     The average time it takes to calculate trapezoidal, Simpson's and scipy.special.dawsn
     with the number of slices returned.
    """

    #intialize variables to store N
    N_trape = 2
    N_simp = 2
    #call the function that returns the difference
    simp = simpson(N_simp,0,4,0)
    trap = trapezoidal(N_trape,0,4,0)
    #calculate the difference between the actual value and the numerical value
    #loop until the difference is of the order given
    while(trap[1] >= (oerror)):
        N_trape *= 2
        trap = trapezoidal(N_trape, 0, 4, 0)
        if (simp[1] >= (oerror)):
            N_simp *= 2
            simp = simpson(N_simp, 0, 4, 0)
    #number of iterations to calculate time average
    n =1000000
    #initialize variables to store total time
    trapeT = 0
    simpT = 0
    specialdwT = 0
    for i  in (1,n):
        t = trapezoidal(N_trape,0,4,1)
        trapeT += t[1]
        s = simpson(N_simp,0,4,1)
        simpT += s[1]
        start1 = time()
        special.dawsn(4)
        end1 = time()
        diff = end1 - start1
        specialdwT += diff
    return np.log2(N_trape), np.log2(N_simp), trapeT/n, simpT/n, specialdwT/n

order9 = ordererror(10e-9)
print("For the Trapezoidal method the number of slices needed for the order of O(-9) is 2^" + str(order9[0]) )
print("For the Simpson method the number of slices needed for the order of O(-9) is 2^" + str(order9[1]))
print("The average time it takes to numerical integrate the trapezoidal method to get an error of 0(-9) is " + str(order9[2]))
print("The average time it takes to numerical integrate the Simpson method to get an error of 0(-9) is " + str(order9[3]))
print("The average time it takes to solve using special.dawsn is "+ str(order9[4]))

#part 2a-iii
#practical estimate of errors
#usingEuler-Maclaurin formula
def practicalerr(N_1):
    """Using the practial estimate of errors from the textbook pg. 153
     --------------------
     :parameter
     N_1: the number of slices, float
     :returns
     An error estimation of the second integral"""
    #Trapezoidal
    trapeI_1 = trapezoidal(N_1, 0, 4, 2)
    trapeI_2 = trapezoidal(N_1*2, 0, 4, 2)
    #Simpson
    simpI_1 = simpson(N_1, 0, 4, 2)
    simpI_2 = simpson(N_1*2, 0, 4, 2)
    return (1/3) * abs(trapeI_2-trapeI_1), (1/15) * abs(simpI_2-simpI_1)
test2 = practicalerr(32)
print("The error estimation for N = 64 for Trapezoidal method is " +str(test2[0]))
print("The error estimation for N = 64 for Simpson's method is " +str(test2[1]))
