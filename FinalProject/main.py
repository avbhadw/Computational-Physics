#final project
import numpy as np
import matplotlib.pyplot as plt
#calculating the ground state of hydrogen using the variational monte carlo method
#then calculate the ground state of heliuem using the same method
#do error analysis stuff.

#hydrogen first
#define a trial function
def trial(r,alpha):
    """computes the trial wavefunction"""
    norm_r = np.linalg.norm(r)
    wave = np.exp(-alpha*norm_r)
    return wave
#a function that calculates the probability density
def prob_density(r, alpha):
    """returns probability density of trial wavefunction"""
    return (trial(r, alpha))**2 #note this is not normalized
#find the local energy
def localE(r,alpha):
    '''compute local energy for the corresponding trial wave function'''
    norm_r = np.linalg.norm(r)
    energy = -1/norm_r-alpha*(alpha-2/norm_r)/2 #where is this mentioned
    return energy

#now use the metropolis algorithm to calculate the ground state of Hydrogen
def Metropolis(N, alpha):
    L = 2
    r = (np.random.rand(3)*2*L-L) #why this method? ?
    E = 0
    E_squared = 0 #for variance of E
    E_avg = 0 #what is this?
    avg = 0 #what is this?
    reject_ratio = 0
    #algorithm
    for i in range(N):
        r_trial = r + 0.2*(np.random.rand(3)*2*L-L) #why
        if prob_density(r_trial, alpha) >= prob_density(r, alpha):
            r = r_trial
        else:
            var = np.random.rand()
            if var < prob_density(r_trial, alpha)/prob_density(r, alpha):
                r = r_trial
            else:
                reject_ratio += 1/N
        E += localE(r, alpha)/N
        E_squared += (localE(r, alpha))**2/N
        E_avg += localE(r, alpha)*-np.linalg.norm(r)/N #what is this
        avg += -np.linalg.norm(r)/N #what is this
    return E, E_squared, E_avg, avg, reject_ratio
#initial parameters
alpha = 0.5 #why
a_iterations = 50
Number_iter = 100 #500
random_walk = 200 #what #200
gamma = 0.5 #what

#initialize arrays to plot graphs
energy_plt = []
alpha_plt = []
variance_plt = []

#define a function to find the minimum values of Energy
for i in range(a_iterations):
    E = 0
    E_square = 0
    dE_dalpha = 0
    Eln = 0
    ln  = 0
    reject_ratio = 0
    for j in range(random_walk): #what is this please explain it
        E_met, E2_met, Eln_met, ln_met, rejections_met = Metropolis(Number_iter, alpha)
        E += E_met/random_walk
        E_square += E2_met/random_walk
        Eln += Eln_met/random_walk
        ln += ln_met/random_walk
        reject_ratio += rejections_met/random_walk  #i have no idea what is being done here
        #define the next alpha
        dE_alpha = 2*(Eln-E*ln)
        alpha = alpha - gamma*dE_dalpha

        #append values to plots
        energy_plt.append(E)
        alpha_plt.append(alpha)
        variance_plt.append(E_square-E**2)
#plotting outside the loop
# plt.plot(energy_plt, alpha_plt)
# plt.show()
#plt.plot(energy_plt)
#plt.show()
#plt.plot(variance_plt)
#plt.show()
fig2 = plt.figure()
ax4 = fig2.add_subplot(111)
plt.title('Hydrogen atom: Energy vs. alpha')
plt.grid()
ax4.plot(alpha_plt, energy_plt, 'ro')
ax4.set_xlabel('Alpha')
ax4.set_ylabel('Energy')
plt.show()

#now trying for helium
