from L09Q02_Functions import *
from numpy import *
from matplotlib import pyplot as plt

############################ CONSTANTS ##########################################
N = 2000
t = 0.01
T = t * N
Lx = 1
Ly = 1
J0 = 1
m = 1
n = 1
c = 1
P = 32
omega = 3.75
Dx = pi * t * c / (2 * Lx)
Dy = pi * t * c / (2 * Ly)
time = arange(0, T + t, t)
# part d
d_omega = 1
omega_0 = arange(0, 9 + d_omega, d_omega)
############################ ARRAYS ##########################################

# coordinate plane
x = arange(0, Lx, P)
y = arange(0, Ly, P)

Hx = zeros((P, P))
Hy = zeros((P, P))
E = zeros((P, P))
Jz = zeros((P, P))

Hx_new = zeros((P, P))
Hy_new = zeros((P, P))
E_new = zeros((P, P))
Jz_new = zeros((P, P))

X_0 = zeros((P, P))
Y_0 = zeros((P, P))
E_hat = zeros((P, P))
J_new = zeros((P, P))

# array to store the values we want to plot
Hx_plt = []
Hy_plt = []
E_plt = []


############################ FUNCTIONS ##########################################
def J_z(w, t):
    J = zeros((P, P))
    for i in range(P):
        for j in range(P):
            x = i * Lx / P
            y = j * Ly / P
            J[i, j] = J0 * (sin(pi * m * x / Lx) * sin(pi * n * y / Ly) * sin(w * t))
    return J


def e_hat(i, j, eh, x, y, j_0):

    var1 = (1 - (j ** 2) * (Dx ** 2) - (i ** 2) * (Dy ** 2)) * eh[i, j]
    var2 = 2 * i * Dy * x[i, j]
    var3 = 2 * j * Dx * y[i, j]
    var4 = j_0[i, j] * 0.01
    var5 = 1 + ((j ** 2) * (Dx ** 2)) + ((i ** 2) * (Dy ** 2))
    var1_4 = var1 + var2 + var3 + var4
    var1_5 = var1_4 / var5
    return var1_5


def fouriercoeff(w):
  
    Hx = zeros((P, P))
    Hy = zeros((P, P))
    E = zeros((P, P))
    for x in time:
        Jz = J_z(w, x)
        J = dst2(Jz)
        X = dHxt2(Hx)
        Y = dHyt2(Hy)
        Ehat = dst2(E)
        # Crank-Nickelson
        for p in range(P):
            for q in range(P):
                E_hat[p, q] = e_hat(p, q, Ehat, X, Y, J)
                X_0[p, q] = X[p, q] - p * Dy * (E_hat[p, q] + Ehat[p, q])
                Y_0[p, q] = Y[p, q] - q * Dx * (E_hat[p, q] + Ehat[p, q])
        # inverse to get the true values
        E = idst2(E_hat)
        Hx = idHxt2(X_0)
        Hy = idHyt2(Y_0)
        # values we want
        Hx_plt.append(Hx[0][P // 2])
        Hy_plt.append(Hy[P // 2][0])
        E_plt.append(E[P // 2][P // 2])
    return Hx_plt, Hy_plt, E_plt

########################## PART C #######################################

hx1, hy1, e1 = fouriercoeff(3.75)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(time, hx1)
plt.title("$H_x$ at (0.5,0)")
plt.ylabel("$H_x$")
plt.xlabel("Time")
plt.show()

plt.plot(time, hy1)
plt.title("$H_y$ at (0,0.5)")
plt.ylabel("$H_y$")
plt.xlabel("Time")
plt.show()

plt.plot(time, e1)
plt.title("$E_z$ at (0.5,0.5)")
plt.xlabel("Time")
plt.ylabel("$E_z$")
plt.show()

# all three plots together
plt.plot(time, E_plt, label='E')
plt.plot(time, Hy_plt, label="$H_y$")
plt.plot(time, Hx_plt, label="$H_x$")
plt.ylabel("Time")
plt.xlabel("$E_z,H_y,H_x$")
plt.title("$E_z,H_y,H_x$")
plt.legend()
plt.show()


########################## PART D #######################################
Emax = 0
e_max = []
for wi in omega_0:
    hx2, hy2, e2 = fouriercoeff(wi)
    e_max.append(max(e2))
    print("1")
#plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(omega_0,e_max)
plt.title("Maximum amplitude of $E_z$ as a function of $\omega$")
plt.xlabel("\omega")
plt.ylabel("$E_z$")
plt.show()

########################## PART E #######################################
omega2e = pi*sqrt(2)
hx0, hy2, e2 = fouriercoeff(omega2e)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(time, hx0)
plt.title("$H_x$ at [0.5,0]")
plt.ylabel("$H_x$")
plt.xlabel("Time")
plt.show()

plt.plot(time, hy2)
plt.title("$H_y$ at [0,0.5]")
plt.ylabel("$H_y$")
plt.xlabel("Time")
plt.show()

plt.plot(time, e2)
plt.title("$E_z$ at [0.5,0.5]")
plt.xlabel("Time")
plt.ylabel("$E_z$")
plt.show()
