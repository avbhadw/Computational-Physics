from numpy import*
import matplotlib.pyplot as plt
#part a
def loaddata():
    data = loadtxt('blur.txt')
    return data
#part a
test = loaddata()
plt.figure(1)
plt.imshow(test,cmap='gray')
plt.title('Blurred image')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('blurredimage.png')
#part b
def GaussSpread(x,y,sig):
    return exp(-(x**2 + y**2)/(2*sig**2))
def GaussSpreadPoints(d):
    rows = d.shape[0]
    cols = d.shape[1]
    valuesgauss = zeros([rows, cols])
    sigma = 25
    for i in range(rows):
        ip = i
        if ip > rows/2:
            ip -= rows #bottom half of rows moved to negative values
        for j in range(cols):
            jp = j
            if jp > cols/2:
                jp -= cols # right half of columns moved to negative
            valuesgauss[i,j] = GaussSpread(ip,jp,sigma) #compute gaussian
    return valuesgauss
#density plot
test = GaussSpreadPoints(test)
plt.figure(2)
plt.imshow(test,cmap='gray')
plt.title('Gaussian Point Spread')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('Gausspointspread.png')
#part c
from numpy.fft import rfft2,irfft2
#call functions to get data
points = loaddata()
gausspoints = GaussSpreadPoints(points)
#calculate the fourier transform
fftpoints = rfft2(points)
fftgauss = rfft2(gausspoints)
#divide one by the other
div = zeros([len(fftgauss), len(fftgauss[0])],dtype=complex)
for i in range(len(fftgauss)):
    for j in range(len(fftgauss[0])):
        if abs(fftgauss[i, j]) > 1e-3:
            div[i, j] = fftpoints[i, j]/fftgauss[i, j]
        else:
            div[i, j] = fftpoints[i, j]
#perform an inverse transfrom to get the unblurred photo
invdiv = irfft2(div)
#plotting
plt.figure(3)
plt.imshow(invdiv, cmap='gray')
plt.title('Deconvoluted Image')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('deconvoluted.png')
