# Question 3 done by Avani
import numpy as np
import matplotlib.pylab as plt
from time import time
def MatrixMultiplication(arr):
    import numpy as np
    import math
    """N is an array that contains the number of iterations
    T stores the time it takes to calculate each matrix multiplication
    A and B are two matrices that will be used for the matrix multiplication 
    C is a matrix that will be used to calculate N^3"""
    t = []#stores n nested loop time
    tdot =[] #stores n np.dot time
    start =np.array([])
    end = np.array([])
    diff = np.array([])
    for n in range(len(arr)):
        A = np.ones([arr[n], arr[n]], float) * 3
        B = np.ones([arr[n], arr[n]], float) * 4
        C =np.zeros([arr[n],arr[n]])
        # stores the matrix multiplication
        x = np.zeros([arr[n], arr[n]])
        # save start time
        start = time()
        # iterating along row of matrix A
        # using direct matrix multiplication
        for i in range(len(A)):
            # iterating along column of matrix B
            for j in range(len(B[0])):
                # iterating along row of matrix B
                for k in range(len(B)):
                    x[i,j] += A[i,k] * B[k,j]
        end = time()
        diff = end-start
        t.append(diff)
        print(t)
        #using npdot
        startdot=time()
        C+=np.dot(A,B)
        enddot =time()
        diffdot =enddot-startdot
        tdot.append(diffdot)
    return t,tdot
iterations = np.arange(0,40,1)
matrix,dot = MatrixMultiplication(iterations)
print("The time is takes to completed a nested matrix multiplication of 200 iterations is "+ str(sum(matrix)))
plt.show()
plt.plot(np.power(iterations,3),matrix)
plt.plot(np.power(iterations,3),dot)
plt.title("Time vs $N^3$")
plt.xlabel("$N^3$")
plt.ylabel("Time")
plt.show()
plt.plot(iterations,matrix)
plt.plot(iterations,dot)
plt.xlabel("N")
plt.ylabel("Time")
plt.title("T vs N")
plt.show()

#log time of N 
# logm = np.log(matrix)
# plt.plot(iterations,logm)
# plt.title("logplot")
#the below program only returns the log plot of time as a function of n^3 for comparison ensurue that the iterations are the same
# def MatrixMultiplicationcube(arr):
#     import numpy as np
#     import math
#     """N is an array that contains the number of iterations
#     T stores the time it takes to calculate each matrix multiplication
#     A and B are two matrices that will be used for the matrix multiplication
#     C is a matrix that will be used to calculate N^3"""
#     t = []#stores n nested loop time
#     cubet = [] #stores cubed time
#     start =np.array([])
#     end = np.array([])
#     diff = np.array([])
#     for n in range(len(arr)):
#         A = np.ones([arr[n], arr[n]], float) * 3
#         B = np.ones([arr[n], arr[n]], float) * 4
#         C=  np.ones([arr[n], arr[n]], float) * 5
#         # stores the matrix multiplication
#         x2 = np.zeros((arr[n],arr[n])) #stores the value of the matrix cubed
#         result =[]
#         # save start time
#         startcube = time()
#         start = time()
#         x = np.zeros([arr[n], arr[n]])
#         # iterating along row of matrix A
#         # using direct matrix multiplication
#         for i in range(len(A)):
#             # iterating along column of matrix B
#             for j in range(len(B[0])):
#                 # iterating along row of matrix B
#                 for k in range(len(B)):
#                     x[i,j] += A[i,k] * B[k,j]
#                 #result.append(x)
#         #print(result)
#         #stuff.append(x)
#         #print(stuff)
#         end = time()
#         diff = end-start
#         t.append(diff)
#         for q in result: #accessing each matrix stored inside the array
#         #N cubed takes the matrix x that is the multiplication of A and B and multiplies it with a new matrix C.
#             for l in range(len(C)):
#         #iterating along the columns of matrix X
#                 for m in range(len(q[0])):
#         #iterating along the row of matrix B
#                     for n in range(len(q)):
#                         x2[l,m] +=C[l,n] * q[n,m]
#         endcube =time()
#         diffcube = endcube-startcube
#         cubet.append(diffcube)
#     return t,cubet
# iterations1 =np.arange(0,200,1)
# matrix1,time1 = MatrixMultiplicationcube(iterations1)
# logt = np.log(time1)
# plt.plot(np.power(iterations1,3),logt)
# plt.title("$N^3$ vs Time")
# plt.xlabel("$N^3$")
# plt.ylabel("Log Time")
# plt.show()


