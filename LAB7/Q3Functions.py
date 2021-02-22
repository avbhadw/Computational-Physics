from numpy import *
from matplotlib.pylab import *
#constants
import scipy.constants as pc
a = pc.physical_constants['Bohr radius'][0]
E0 = pc.physical_constants['Rydberg constant times hc in eV'][0]
m = pc.m_e# electron mass
hbar = pc.hbar
e = pc.e# elementary charge
eps0 = pc.epsilon_0
# stepsize and r-inf
h = 0.0002 * a
rinf = 20 * a
r = arange(h, rinf, h)

#normalizing the wave function
def Simpson(r_a):
    N = len(r_a)
    r = abs(r_a)
    sum = r[-1]+r[0]
    for i in range(1, N, 2):
        sum += 4*r[i]**2
    for i in range(2, N, 2):
        sum += 2*r[i]**2
    return (h/3*a)*sum
# Potential function
def V(x):
    vr = -(e**2) / (4 * pi * eps0 * x)
    return vr
#writing the given equation in 2 first order ODEs
def f(r, x, E, l):
    R = r[0]
    S = r[1]
    fR = S / x ** 2
    fS = ((2 * m * x ** 2 / hbar ** 2) * (V(x) - E) + l * (l + 1)) * R
    return array([fR, fS], float)


# Calculate the wave function for a given energy
def solve(E, l):
    R = 0.0
    S = 1.0
    global wavefunc
    wavefunc = zeros((len(r), 2), float)
    #initial values
    wavefunc[0] = array([R, S])

    for i in arange(len(wavefunc) - 1):
        k1 = h * f(wavefunc[i], r[i], E, l)
        k2 = h * f(wavefunc[i] + 0.5 * k1, r[i] + 0.5 * h, E, l)
        k3 = h * f(wavefunc[i] + 0.5 * k2, r[i] + 0.5 * h, E, l)
        k4 = h * f(wavefunc[i] + k3, r[i] + h, E, l)
        wavefunc[i + 1] = wavefunc[i] + (k1 + 2 * k2 + 2 * k3 + k4)/6
    return wavefunc
# Main program to find the energy using the secant method
def secant(n,l):
    E1 = (-15 * e / n ** 2)
    E2 = (-13 * e / n ** 2)
    R2 = solve(E1, l)[-1, 0]

    target = e / 1000
    while abs(E1 - E2) > target:
        R1, R2 = R2, solve(E2, l)[-1, 0]
        E1, E2 = E2, E2 - R2 * (E2 - E1) / (R2 - R1)
    return E2
#analytic solutions
def analytic(type):
    if (str(type) == str('R1')):
        num = 1*exp(-r/a)
        dem = sqrt(pi)*(a ** (3/2))
        r1 = num/dem
        #nomalize the function
        norm_r1 = sqrt(Simpson(r1)*(a**(1.5)))
        return r1/norm_r1
    elif(str(type)== 'R20'):
        x = 1/(4*sqrt(2*pi)*(a**1.5))
        y = (2-(r/a))*exp(-r/(2*a))
        r20 = x*y
        #normalize
        norm_r20 = sqrt(Simpson(r20)*(a**(1.5)))
        return r20/norm_r20
    else:
        x = 1 / (4 * sqrt(2 * pi) * (a ** 2.5))
        y = (r/a)*exp(-r/(2*a))
        r21 = x*y
        #normalize
        norm_r21 = sqrt(Simpson(r21)*(a ** 1.5))
        return r21/norm_r21

