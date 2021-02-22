from numpy import *
from matplotlib.pyplot import *
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
def shift(rx, ry):
    rx_shifted = zeros(len(rx))
    ry_shifted = zeros(len(rx))
    for a in range(-1, 2):
        rx_shifted = (rx + 4 * a)
        for b in range(-1, 2):
            ry_shifted = (ry + 4 * b)
    return rx_shifted, ry_shifted

def VerletAlgo(N, h, stepno):
    """Update the x and y positions using the Verlet method """
    # initial arrays
    # velocity
    global initvy, initvx
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
    #initialize arrays with given values
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
                vx2 += 1 / 2 * h * a(x[0, j] - x[0, k], y[0, j] - y[0, k], "x")
                vy2 += 1 / 2 * h * a(x[0, j] - x[0, k], y[0, j] - y[0, k], "y")
        # moving columnwise for updating
        return vx2, vy2
    def k_update(x,y):
        k1 = zeros(N)
        k2 = zeros(N)
        for m in range(0, N):
            for n in range(0, N):
                k1[m] += h * a(x[i + 1, m] - x[i + 1, n], y[i + 1, m] - y[i + 1, n], 'x')
                k2[m] += h * a(x[i + 1, m] - x[i + 1, n], y[i + 1, m] - y[i + 1, n], 'y')
        return k1, k2

    for p in range(0, 9):
        rx, ry = shift(rx, ry)
        initvx, initvy = vhalf(rx, ry, vx, vy, h)
        for i in range(0, stepno - 1):
            # implementing the boundary and periodic condition
            rx[i + 1] = rx[i] + (initvx * h)
            rx[i + 1] = mod(rx[i + 1], Lx)
            ry[i + 1] = ry[i] + (initvy * h)
            ry[i + 1] = mod(ry[i + 1], Ly)
        # update k, eq9
            kx, ky = k_update(rx, ry)
        # update velocity
            vx[i + 1] = vx[i] + (0.5 * kx)
            vy[i + 1] = vy[i] + (0.5 * ky)
        # update the last equation 10
            initvx += kx
            initvy += ky
            time[i + 1] = time[i] + h
        return rx, ry, time
x, y, t = VerletAlgo(16,0.01, 1000)
# plotting the trajectory
figure()
rc('text', usetex=True)
rc('font', family='serif')
plot(x, y, marker='.')
title('Periodic and Boundary Conditions ')
xlabel('$x$', fontsize=12)
ylabel('$y$', fontsize=12)
show()

