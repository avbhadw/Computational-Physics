#Question 2B
#exercise 5.4 in the textbook
#part a: Numerical Integration for the Bessel Function
import numpy as
def Bessel(m, x):
    def f(m, x, ang):
        return 1/np.pi*(np.cos(m * ang - x * np.sin(ang)))
    #integrating using simpson rule
    #start and end point
    a = 0.0
    b = math.pi
    #steps
    N = 1000
    h = (b - a) / N
    odd = 0
    for k in range(1, N, 2):
        odd += 4.0 * f(m, x, a + k * h)
    even = 0
    for k in range(2, N, 2):
        even += 2.0 * f(m, x, a + k * h)
    return (h / 3.0) * (even + odd + f(m, x, a) + f(m, x, b))

#2b-i

x =np.arange(0,21,0.1)
#initialize array to store data
j0 = np.zeros(len(x))
j1 = np.zeros(len(x))
j2= np.zeros(len(x))
for i in range(len(x)):
    j0[i] = (Bessel(0, x[i]))
    j1[i] = (Bessel(1, x[i]))
    j2[i] = (Bessel(2, x[i]))

#plotting
#storing values from special.jv
jv0 = special.jv(0, x)
jv1 = special.jv(1, x)
jv2 = special.jv(2, x)
#storing difference between true value and numerical
diff0 = abs(j0 - jv0)
diff1 = abs(j1 - jv1)
diff2 = abs(j2 - jv2)
#plot numerical vs x
plt.plot(x,jv0,label = "$J_0$")
plt.plot(x,jv1,label = "$J_1$")
plt.plot(x,jv2,label ="$J_2$")
plt.xlabel("X from 0 to 20")
plt.ylabel("scipy.special.jv")
plt.title("X as a function of special.jv")
plt.legend()
plt.show()
#plot special.jv vs x
plt.plot(x, j0, label="$J_0$")
plt.plot(x, j1, label="$J_1$")
plt.plot(x, j2, label="$J_2$")
plt.xlabel("X from 0 to 20 ")
plt.ylabel("Numerical Bessel function ")
plt.title("X as a function of Numerical Bessel function")
plt.legend()
plt.show()
#plot diff between special.jv and numerical integration
plt.plot(x,diff0,label ="Difference $J_0$")
plt.plot(x,diff1,label ="Difference $J_1$")
plt.plot(x,diff2,label = "Difference $J_2$")
plt.suptitle("Difference between special.jv and numerical integration ", y = 0.98)
plt.xlabel("X from 0 to 20 ")
plt.ylabel("Difference")
plt.legend()
plt.show()

#part b
def I(radius):
    """ density plot of the intensity of the diffraction pattern  """
    intensities = np.zeros(([len(x),len(x)]))
    for i in range (len(radius)):
      for j in range (len(radius)):
          if (radius[i,j]) == 0.000001:
              intensities[i,j] = 0.25
          else:
              Lambda = 0.5  # in micrometers
              k = 2 * np.pi / Lambda * radius[i, j]
              intensities[i,j] =  (Bessel(1, k) / k) ** 2
    print(intensities)
    plt.pcolormesh(x,y,intensities, vmax=0.001)
    plt.suptitle("Intensity of light $W/m^2$")
    plt.xlabel("x ($\mu m$)")
    plt.ylabel("y ($\mu m$)")
    plt.colorbar()
    plt.gray()
    plt.show()
    return intensities
x1 = np.arange(-1, 1, 0.01)
y1 = np.arange(-1, 1, 0.01)
x, y = np.meshgrid(x1, y1)
radius = np.sqrt(x ** 2 + y ** 2)
intensities = I(radius)

