from numpy import*
from matplotlib.pyplot import *
#now im working on q3
#2b
#define acceleration as components
def a(x,y,type):
    r = x**2+y**2
    a_x = 12 * ((2 * x / (r ** 7)) - (x / (r ** 4)))
    a_y = 12 * ((2 * y / (r ** 7)) - (y / (r ** 4)))
    if (str(type) == 'x' or 'y') and r == 0:
        return 0
    if str(type) == 'x':
        return a_x
    if str(type) == 'y':
        return a_y
    else:
        return a_x,a_y
# def potentialE(x,y):
#     r = x**2+y**2
#     pe = 0
#     for i in range(len(x)):
#         for j in range(len(y)):
#             if i ! = j:
#                 pe +=4*1*(1/r**12 - 1/r**6)
#         return pe
def kinecticE(v_1,v_2):
    return 0.5*v_1**2+0.5*v_2**2

def VerletAlgo(N,start,h,stepno):
    """Update the x and y positions using the Verlet method """
    #initial arrays
    #velocity
    global initvy, initvx
    vx = zeros([stepno, N])
    vy = zeros([stepno, N])
    #position
    rx = zeros([stepno, N])
    ry = zeros([stepno, N])
    #time
    time = zeros(stepno)
    time[0] = start
    Lx = 4.0
    Ly = 4.0
    dx = Lx/sqrt(N)
    dy = Ly/sqrt(N)
    x_grid = arange(dx/2,Lx,dx)
    y_grid = arange(dy/2,Ly,dy)
    xx_grid, yy_grid = meshgrid(x_grid,y_grid)
    x_initial = xx_grid.flatten()
    y_initial = yy_grid.flatten()
    for i in range(N):
        rx[0][i] = x_initial[i]
        ry[0][i] = y_initial[i]
        vx[0][i] = 0
        vy[0][i] = 0
    #calculate the distance rij to particle j from i for all j
    #print(rx)
    # for i in range(N):
    #     for j in range(N):
    #         if i != j:
    #             distx = x_initial[i] - x_initial[j]
    #             disty = y_initial[i] - y_initial[j]
    #             #calculating accleration at each point
    #             acc_x,acc_y = a(distx,disty,'l')
    #             #print(acc_x)
    # #print(dist)

    #implementing verlet algorithm
    def vhalf(x,y,v_x,v_y,h):
        vx2 = v_x[0]
        vy2 = v_y[0]
        for j in range(0, N):
            for k in range(0, N):
                vx2 += 1 / 2 * h * a(x[0, j] - x[0, k], ry[0, j] - y[0, k], "x")
                vy2 += 1 / 2 * h * a(x[0, j] - x[0, k], ry[0, j] - y[0, k], "y")
        return vx2, vy2
    for i in range(0, stepno-1):
        initvx,initvy = vhalf(rx, ry, vx, vy, h)
        #update position
        rx[i+1] = rx[i] + (initvx * h)
        ry[i+1] = ry[i] + (initvy * h)
        # update k
        kx = zeros(N)
        ky = zeros(N)
        for m in range(0, N):
            for n in range(0, N):
                kx[m] += h * a(rx[i + 1, m] - rx[i+1, n], ry[i + 1, m] - ry[i+1, n], 'x')
                ky[m] += h * a(rx[i + 1, m] - rx[i+1, n], ry[i + 1, m] - ry[i+1, n], 'y')
        #update velocity
        vx[i+1] = vx[i]+(0.5*kx)
        vy[i+1] = vy[i]+(0.5*ky)
        #update the last equation
        initvx += kx
        initvy += ky
        #computing potential and kinectic energy
        potential = 4*1*(1/(rx[i]**2+ry[i]**2)**12 - 1/(rx[i]**2+ry[i]**2)**6)
        #computing the kinetic energy
        kinectic = 0.5*vx[i]**2+0.5*vy[i]**2
        #update time
        time[i + 1] = time[i] + h
    return rx, ry, time,potential,kinectic
#x,y,t = VerletAlgo(2,[[4.5,4],[5.2,4]],[[0,0],[0,0]],0,0.01,100)
x,y,t,pe,ke = VerletAlgo(16,0,0.01,1000)
#plotting the trajectory
#print(x)
plot(x,y,marker='.')
#gca().set_aspect('equal', adjustable='box')
show()
plot(pe,t)
plot(ke,t)
plot(pe+ke,t)
show()
#plot(t,x)
#show()
#rough stuff just plotting
#now for plotting running a loop for three times
# #r1
# r1 = [[[4, 4],[5.2,4]], [[4.5, 4],[5.2,4]], [[2, 3],[3.5,4.4]]]
# r2 = [[5.2, 4], [5.2, 4], [3.5, 4.4]]
# initialv = [[0, 0], [0, 0]]
# initialt = 0
# stepsize = 0.01
# numberStep = 100
# plots = zeros(3)
# for i in range(0, 2):
#     plots[i] = VerletAlgo(r1[i], initialv, initialt, stepsize, numberStep)
#     plot(plots[0], plots[1], marker='.')
#     savefig(str(i)+'.png')
#     show()