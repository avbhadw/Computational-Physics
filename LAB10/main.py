from numpy import*
from matplotlib.pyplot import *

#3a
a = 0
b = 1
N = 10000
n = 100
I_imp = zeros(n)
I_mean = zeros(n)
#using Importance sampling #change this into a function and call this function in the loop
def f(x):
    return 1/(sqrt(x)*(1+exp(x)))
def weight(x):
    return (1/sqrt(x))
for i in range(0, 100):
    x1 = random.rand(N)**2
    x2 = random.rand(N)
    I_imp[i] = sum(f(x1)/weight(x1))
    I_mean[i] = sum(f(x2))
#plotting
I_imp = 2*I_imp/N
rc('text', usetex=True)
rc('font', family='serif',size=12)
hist(I_imp, 10, range=[0.8, 0.88])
title('Importance Sampling')
xlabel('Result')
ylabel('Counts')
show()
I_mean = I_mean/N
rc('text', usetex=True)
rc('font', family='serif',size=12)
hist(I_mean, 10, range=[0.8, 0.88])
title('Mean Value Method')
xlabel('Results')
ylabel('Counts')
show()

#3b
def g(x):
    return exp(-2*abs(x-5))
def weighted(x):
    return exp(-(x-5)**2/2)*1/sqrt(2*pi)
#start and end of the integral
a = 0
b = 10

#initalize arrays to store values
I_mean = zeros(n)
I_imp = zeros(n)

#calculate the values
for i in range(0, 100):
    x1 = random.normal(5, 1, N)
    I_imp[i] = sum(g(x1)/weighted(x1))
    x2 = random.rand(N)*(b-a)
    I_mean[i] = sum(g(x2))

#plotting
I_imp = I_imp/N
rc('text', usetex=True)
rc('font', family='serif',size=12)
hist(I_imp, 10, range=[0.95, 1.05])
title('Importance Sampling')
xlabel('Result')
ylabel('Counts')
show()
I_mean = (b-a)*I_mean/N
rc('text', usetex=True)
rc('font', family='serif',size=12)
hist(I_mean, 10, range=[0.95, 1.05])
title('Mean Value Method')
xlabel('Results')
ylabel('Counts')
show()

