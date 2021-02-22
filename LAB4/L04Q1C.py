from SolveLinear import PartialPivot
from numpy import*
import matplotlib.pyplot  as plt
# define constants
C1 = 1*10**-6
C2 = 0.5*10**-6
x_p = 3
w = 1000 #s^-1

#define function that declares A and v and returns x
def matrix(R6):
    R1 = R3 = R5 = 1000  #  ohms
    R2 = R4 = 2000  #  ohms
    A = array([[1 / R1 + 1 / R4 + 1j * w * C1,  -1j * w * C1,  0],
        [-1j * w * C1,  1 / R2 + 1 / R5 + 1j * w * C1 + 1j * w * C2, -1j * w * C2],
        [0, -1j * w * C2,  1 / R3 + 1 / R6 + 1j * w * C2]], complex)
    v = array([x_p / R1, x_p / R2, x_p / R3], complex)
    return PartialPivot(A,v)
#R6 = 2000
x = matrix(2000)
#for plotting
period = linspace(0, 4 * pi / w, 101)
#finding the values V
for i in range(len(x)):
    print('The magnitude of V'+str(i)+' is:' + str(abs(x[i]))+',The phase of V'+str(i)+' is:' +str(180/pi*angle(x[i])))
    #plotting
    plt.plot(period, real(x[i]*exp(1j*w*period)), label='$V_'+str(i+1)+'$')
plt.legend()
plt.title('$R_6 = 2000\Omega$')
plt.xlabel('Time(s)')
plt.ylabel('Voltage (V)')
plt.show()
plt.savefig('R6')
# using an inductor
xL = matrix(2000j)
for i in range(len(xL)):
    print('The magnitude of V' + str(i) + ' is:' + str(abs(xL[i]))+',The phase of V' + str(i) + ' is:' + str(180/pi*angle(xL[i])))
    plt.plot(period,real(xL[i]*exp(1j*w*period)),label = '$V_'+str(i+1)+'$')
plt.legend()
plt.title('Replacing $R_6$ with an Inductor')
plt.xlabel('Time (s)')
plt.ylabel('Voltage(V)')
plt.show()
plt.savefig('Inductor_iR6')