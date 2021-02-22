import numpy as np
import math
from scipy import special
from gaussxw import gaussxw
def f(t):
    "Returns the dawson function"
    return math.exp(-(4**2))*(math.exp(t ** 2))
def f1(u_10):
    T_a = np.arange(-40,15,1)
    delta = 4.3+0.145*T_a+0.00196*(T_a**2)
    t_h = np.array([24, 48, 72])
    mean_u = 11.2+0.365*T_a+0.00706*(np.power(T_a,2))+0.9*np.lo(t_h)
    return (1/np.sqrt(2*np.pi)*delta)*((np.exp(mean_u-u_10)**2)/2*delta**2)
def trapezoidal(N,a,b):
    """This function calculates the trapezoidal method
    --------------------
    :parameter
    N: the number of slices, float
    a : Start value, float
    b: End value, float
    T: which values to return in the function, boolean
    --------------------"""

    #calculate the width between two points
    h=(b-a)/N
    #store values for f(a) and f(b)
    s = 0.5 * f(a) + 0.5 * f(b)


    for k in range(1, N):
        s += f(a + k * h)
    sol = h*s


    #difference between true value and numerical approximation
    compareT= abs(sol - special.dawsn(b))
    return sol

def simpson(N,a,b):
    """This function calculates the trapezoidal method
       --------------------
       :parameter
       N: the number of slices, float
       a : Start value, float
       b: End value, float
       T: which values to return in the function, boolean
       --------------------"""

    #width of each slice
    h = (b - a) / N

    even =0
    odd = 0

    #for odd values
    for k in range(1,N,2):
        odd += 4 * f(a+k*h)
    #for even values
    for k in range (2,N,2):
        even += 2 * f(a+k*h)
    integral = (1/3*h*(f(a)+f(b)+odd+even))
    #stores the time
    #stores the difference between the true value and numerical approximation
    compareS= abs(integral - special.dawsn(b))
    return integral
def gaussian(N,a,b):
    # Calculate the sample points and weights, then map them
    # to the required integration domain
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    # Perform the integration
    s = 0.0
    for k in range(N):
        s += wp[k]*f1(xp[k])

    return s