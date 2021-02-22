#q1 comparing the three methods
import numpy as np
import math
from scipy import special
import matplotlib.pyplot  as  plt
from L03Functions import trapezoidal,simpson,gaussian

N = np.arange(8,2048,1)
trape = []
simp = []
gauss = []
errorgauss = []
errorpracgauss = []
errortrape =[]
errorpractrape=[]
errorsimp = []
errorpracsimp = []
# for n in range(8,2048):
#     trape.append(trapezoidal(n,0,4))
#     simp.append(simpson(n,0,4))
#     gauss.append(gaussian(n,0,4))
#     #calcuating the relative error compared with the true val
#     errorgauss.append(abs(gauss[-1] - special.dawsn(4)))
#     errortrape.append(abs(trape[-1] - special.dawsn(4)))
#     errorsimp.append(abs(simp[-1] - special.dawsn(4)))
#     #calculate the relative error using practical estimate
#     errorpracgauss.append((gaussian(2*n,0,4) - gaussian(n,0,4)))
#     errorpractrape.append(1/3 * (gaussian(2*n,0,4) - gaussian(n,0,4)))
#     errorpracsimp.append(1/15 * (gaussian(2*n,0,4) - gaussian(n,0,4)))
# #plotting the magnitude of relative error against N
# plt.plot(N, errorgauss, label="Gauss")
# plt.plot(N, errorsimp, label="Simpson")
# plt.plot(N, errortrape, label="Trapezoidal")
# plt.legend()
# plt.title("Relative error compared with true value")
# plt.xlabel("N")
# plt.ylabel("Error")
# #plt.loglog()
# plt.show()
# plt.plot(N, errorpracgauss, label="Gauss")
# plt.plot(N, errorpracsimp, label="Simpson")
# plt.plot(N, errorpractrape, label="Trapezoidal")
# plt.loglog()
# plt.legend()
# plt.title("Practical estimate of errors")
# plt.xlabel("N")
# plt.ylabel("Error")
# plt.show()

#2b-i
def blowingsnow(u_10,t_h):
    T_a = np.arange(-40,15,1)
    delta = 4.3+0.145*T_a+0.00196*(T_a**2)
    mean_u = 11.2+0.365*T_a+0.00706*(T_a**2)+0.9*np.log(i)
    return (1/np.sqrt(2*np.pi)*delta)*((np.exp(mean_u-u_10)**2)/2*delta**2)
u_10 = np.array([6,8,10])
t_h = np.array([24, 48, 72])
T_a = np.arange(-40, 15, 1)
probability = np.zeros((24,25))
snow = []
for i in range(len(u_10)):
    print(i)
    for j in t_h:
        probability[i,j] = blowingsnow(u_10[i],j)
        print(probability[i,j])
