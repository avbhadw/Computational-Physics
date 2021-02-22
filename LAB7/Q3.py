from numpy import*
from matplotlib.pyplot import *
from Q3Functions import*
import scipy.constants as pc
e = pc.e# elementary charge
E0 = pc.physical_constants['Rydberg constant times hc in eV'][0]
#part b
print('The numerical value is', secant(1, 0)/e)
print('The theorectical value is, ', -E0)
print('The numerical value is', secant(2, 0)/e)
print('The theoretical value is', -E0/4)
print('The numerical value is', secant(2, 1)/e)
print('The theoretical value is', -E0/4)
#n = 1 l = 0
figure(1)
r_1 = solve(secant(1, 0), 0)[:, 0]
norm = sqrt(Simpson(r_1)*(a**(1.5)))
rc('text', usetex=True)
rc('font', family='serif')
plot(r/a, r_1/norm, label="numerical", linestyle='dashed')
plot(r/a, analytic('R1'), label='analytic', linestyle='dashed')
title('Radial Function for n = 1, $l$ = 0')
ylabel('Wavefunction')
xlabel('Bohr Radii ($5.2918\\times 10^{-11} $)')
legend()
savefig('n1l0.png')
#n = 2 l = 0
figure(2)
r_2 = solve(secant(2, 0), 0)[:, 0]
norm2 =sqrt(Simpson(r_2)*(a**(1.5)))
rc('text', usetex=True)
rc('font', family='serif')
plot(r/a, r_2/norm2, label="numerical", linestyle='dashed')
plot(r/a, analytic('R20'), label='analytic', linestyle='dashed')
title('Radial Function for n = 2, $l$ = 0')
ylabel('Wavefunction')
xlabel('Bohr Radii ($5.2918\\times 10^{-11} $)')
axhline(y=0, xmin=0.05, xmax=1, color='red', label='x = 0 ')
legend()
savefig('n2l0.png')
# n  2 l = 1
figure(3)
r_21 = solve(secant(2, 1), 1)[:, 0]
norm21 =sqrt(Simpson(r_21)*(a**(1.5)))
rc('text', usetex=True)
rc('font', family='serif')
plot(r/a, r_21/norm21, label="numerical", linestyle='dashed')
plot(r/a, analytic('R21'), label='analytic', linestyle='dashed')
title('Radial Function for n = 2, $l$ = 0')
ylabel('Wavefunction')
xlabel('Bohr Radii ($5.2918\\times 10^{-11} $)')
#axhline(y=0, xmin=0.05, xmax=1, color='red', label='x = 0 ')
legend()
savefig('n2l1.png')