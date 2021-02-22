import numpy as np
import matplotlib.pyplot as plt

def trial_wf(r1,r2,alpha):
    norm_r1 = np.linalg.norm(r1)
    norm_r2 = np.linalg.norm(r2)
    r_12 = np.linalg.norm(r1-r2)
    wave_func = np.exp(-2*norm_r1)*np.exp(-2*norm_r2)*np.exp(r_12/(2*(1+alpha*r_12)))
    return wave_func
def probability_den(r1,r2,alpha):
    return trial_wf(r1,r2,alpha)**2
def local_E(r1,r2,alpha):
    norm_r1 = np.linalg.norm(r1)
    norm_r2 = np.linalg.norm(r2)
    r_12 = np.linalg.norm(r1 - r2)
    dot = np.dot(r1/norm_r1 - r2/norm_r2, r1-r2)
    energy = -4+dot/(r_12*(1+alpha*r_12)**2)-1/(r_12*(1+alpha*r_12)**3)-1/(4*(1+alpha*r_12)**4)+1/r_12
    return energy
def Metropolis(N,alpha):
    L = 1
    r1 = np.random.rand(3)*2*L-L #3 degrees of freedom
    r2 = np.random.rand(3)*2*L-L
    E = 0

    step = 0
    max_step = 500
    for i in range(N):
        num = np.random.rand() #choose a random number
        step += 1
        if num < 0.5:
            r1_trial = r1+0.5*(np.random.rand(3)*2*L-L)
            r2_trial = r2
        else:
            r2_trial = r2+0.5*(np.random.rand(3)*2*L-L)
            r1_trial = r1
        if probability_den(r1_trial,r2_trial,alpha) >= probability_den(r1,r2,alpha):
            r1 = r1_trial
            r2 = r2_trial
        else:
            num2 = np.random.rand()
            var = probability_den(r1_trial,r2_trial,alpha)/probability_den(r1,r2,alpha)
            if num2<var:
                r1 = r1_trial
                r2 = r2_trial

        if step > max_step:
            div = N-max_step
            localE = local_E(r1,r2,alpha)
            E += localE/(div)
            r_12 = np.linalg.norm(r1-r2)

    return E
#initialize values
alpha = 0
alpha_iterations = 6
N_met = 5000
ramdom_walker = 200

energy_plt = []
alpha_plt = []

#minimize energy and plot

for i in range(alpha_iterations):
    E = 0
    for j in range(ramdom_walker):
        E_met = Metropolis(N_met, alpha)
        E += E_met /ramdom_walker
    alpha = alpha + 0.05
    #append plotting arrays
    energy_plt.append(E)
    alpha_plt.append(alpha)

print(min(energy_plt))
fig2 = plt.figure()
ax4 = fig2.add_subplot(111)
plt.title('Helium atom: Energy vs. alpha')
plt.grid()
ax4.plot(alpha_plt, energy_plt, 'ro')
ax4.set_xlabel('Alpha')
ax4.set_ylabel('Energy')
plt.show()