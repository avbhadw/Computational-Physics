# This program calculates the total energy and magnetization
# for a 1D Ising model with N dipoles
# Author: Nico Grisouard, University of Toronto
# Date: 24 November 2020

# import modules
import numpy as np
from random import random, randrange
from matplotlib import pyplot as plt
import matplotlib.animation as animation
#still need to figure this out


def energyfunction(J_, dipoles):
    """ function to calculate energy """
    row = np.sum(dipoles[0:-1, :]*dipoles[1:, :])
    col = np.sum(dipoles[:, 0:-1]*dipoles[:, 1:])
    return -J_*(row+col)


def acceptance(Enew, E, T):
    """ Function for acceptance probability; to be completed """
    # Do stuff here
    return True if random() < np.exp((E-Enew)/(T*kB)) else False


# define constants
kB = 1.0
T = 3.0
J = 1.0
num_dipoles = 100
l = 20
N = 100000

# generate array of dipoles and initialize diagnostic quantities
dipoles = np.ones([l, l], int)  # hint: this will not work
for i in range(l):
    for j in range(l):
        dipoles[i, j] = randrange(0, 2)
        if dipoles[i, j] == 0:
            dipoles[i, j] = -1

energy = []  # empty list; to add to it, use energy.append(value)
magnet = []  # empty list; to add to it, use magnet.append(value)

E = energyfunction(J, dipoles)
energy.append(E)
magnet.append(np.sum(dipoles))

#Empty plot for annimation stuff
framerate = 1000
fig = plt.figure(0)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Evolution of system with T='+str(T))
plot = []
for i in range(N):
    x_picked = randrange(l)  # choose a victim
    y_picked = randrange(l)
    dipoles[x_picked, y_picked] *= -1  # propose to flip the victim
    Enew = energyfunction(J, dipoles)  # compute Energy of proposed new state
    # store energy and magnetization

    if acceptance(Enew, E, T):
        E = Enew
    else:
        dipoles[x_picked, y_picked] *= - 1

    if i % framerate == 0:
        plot.append([plt.imshow(dipoles, cmap='binary', animated='True')])
    energy.append(E)
    magnet.append(np.sum(dipoles))
#annimation stuff
anni = animation.ArtistAnimation(fig, plot, interval=2, blit=True)
FFwriter = animation.FFMpegWriter(fps=30)
anni.save('test.mp4', writer=FFwriter)


# plot energy, magnetization

#energy plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=12)
plt.figure()
plt.plot(energy)
plt.title('Total Energy for T = 3')
plt.xlabel('Monte Carlo steps ')
plt.ylabel('Energy')
plt.show()

#magnet plot
plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=12)
plt.plot(magnet)
plt.title("Total Magnetization for T = 3")
plt.xlabel("Monte Carlo steps")
plt.ylabel("Magnetization")
plt.show()


