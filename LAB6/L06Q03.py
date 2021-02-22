from numpy import *
from matplotlib.pyplot import *
# now im working on q3
# define acceleration as components
def a(x, y, type):
    # calculate distance r
    r = x ** 2 + y ** 2
    # calculate components
    a_x = 12 * ((2 * x / (r ** 7)) - (x / (r ** 4)))
    a_y = 12 * ((2 * y / (r ** 7)) - (y / (r ** 4)))
    if (str(type) == 'x' or 'y') and r == 0:
        return 0
    if str(type) == 'x':
        return a_x
    if str(type) == 'y':
        return a_y
    else:
        return a_x, a_y

def PotentialE(x,y):
    PEx = 0
    for m in range(0, len(x)):
        for n in range(0, len(y)):
            if m != n:
                if m > n:
                    distx = x[n] - x[m]
                    disty = y[n] - y[m]
                    dist = (distx ** 2 + disty ** 2)
                    PEx += 2 * ((1 / (dist ** 6)) - (1 / (dist ** 3)))
    return PEx
def KinecticE(vel_x, vel_y):
    KE = 0
    for i in range(len(vel_x)):
        KE += vel_x[i]**2
        KE += vel_y[i]**2
    return KE
def VerletAlgo(N, h, stepno):
    """Update the x and y positions using the Verlet method """
    # initial arrays
    # velocity
    global initvy, initvx, PE, KE
    vx = zeros([stepno, N])
    vy = zeros([stepno, N])
    # position
    rx = zeros([stepno, N])
    ry = zeros([stepno, N])
    # time
    time = zeros(stepno)
    Lx = 4.0
    Ly = 4.0
    dx = Lx / sqrt(N)
    dy = Ly / sqrt(N)
    x_grid = arange(dx / 2, Lx, dx)
    y_grid = arange(dy / 2, Ly, dy)
    xx_grid, yy_grid = meshgrid(x_grid, y_grid)
    x_initial = xx_grid.flatten()
    y_initial = yy_grid.flatten()

    PotE = zeros(stepno)
    KinE = zeros(stepno)
    # initialize the variables with given values
    for i in range(N):
        rx[0][i] = x_initial[i]
        ry[0][i] = y_initial[i]
        vx[0][i] = 0
        vy[0][i] = 0
    # implementing verlet algorithm
    # this is the eq 7
    def vhalf(x, y, v_x, v_y, h):
        vx2 = v_x[0]
        vy2 = v_y[0]
        for j in range(0, N):
            for k in range(0, N):
                vx2 += 1 / 2 * h * a(x[0, j] - x[0, k], ry[0, j] - y[0, k], "x")
                vy2 += 1 / 2 * h * a(x[0, j] - x[0, k], ry[0, j] - y[0, k], "y")
        # moving columnwise for updating
        return vx2, vy2
    def kupdate(x,y,h,it):
        kx = zeros(N)
        ky = zeros(N)
        i = it
        for m in range(0, N):
            for n in range(0, N):
                kx[m] += h * a(x[i + 1, m] - x[i + 1, n], y[i + 1, m] - y[i + 1, n], 'x')
                ky[m] += h * a(x[i + 1, m] - x[i + 1, n], y[i + 1, m] - y[i + 1, n], 'y')
        return kx,ky
    for i in range(0, stepno - 1):
        initvx, initvy = vhalf(rx, ry, vx, vy, h)
        # update position, eq8
        rx[i + 1] = rx[i] + (initvx * h)
        ry[i + 1] = ry[i] + (initvy * h)
        # update k, eq9
        kx, ky = kupdate(rx, ry, h, i)
        # update velocity
        vx[i + 1] = vx[i] + (0.5 * kx)
        vy[i + 1] = vy[i] + (0.5 * ky)
        # update the last equation 10
        initvx += kx
        initvy += ky
        #call functions for KE and PE
        # computing potential and kinectic energy
        PotE[i] = PotentialE(rx[i], ry[i])
        # computing the kinetic energy
        KinE[i] = KinecticE(vx[i], vy[i])
        # update time
        time[i + 1] = time[i] + h
    return rx, ry, time, KinE, PotE


x, y, t, KE, PE = VerletAlgo(16, 0.01, 1000)

E = KE + PE
E_diff = zeros(999)
for i in range(1, 1000):
    E_diff[i - 1] = 100 * (abs((E[i] - E[i - 1]) / E[i]))
# plotting the trajectory
# print(x)
figure()
plot(x, y, marker='.')
gca().set_aspect('equal', adjustable='box')
show()

figure()
plot(t, KE, label='KE')
plot(t, PE, label='PE')
plot(t, E, label='Total energy')
legend()
show()

figure()
scatter(t[1:], E_diff)
title('% Energy change from each step vs time')
show()
