from numpy import*
import matplotlib.pyplot as plt
from dft import dft
from scipy.signal import find_peaks
#read sunspots.txt
x, y= loadtxt("sunspots.txt",unpack='True')
y -= mean(y)
plt.figure(1)
plt.plot(x,y)
plt.title('Sunspots as a function of time')
plt.xlabel('Months')
plt.ylabel('Sunspot number')
plt.savefig('sunspots.png')
#plotting c_k^2
fourier = dft(y)
fftsq = abs(fourier**2)
plt.figure(2)
plt.plot(fftsq)
plt.title('$| c_k ^2|$ as a function of k')
plt.xlabel('k')
plt.ylabel('$| c_k^2|$')
plt.savefig('fouriersunspots.png')
#zooming into  plot to find k
plt.figure(3)
plt.plot(fftsq)
plt.xlim(0,50)
plt.title('$| c_k^2 |$ as a function of k ')
plt.ylabel('$|c_k^2|$')
plt.xlabel('k')
plt.savefig('zoomfft.png')
plt.show()

