#Question 1 done by Avani
import numpy as np
import matplotlib.pylab as plt
import math
def orbit(initx,inity,initvx,initvy,step,type,name):
    """This is a function that will plot the orbit of a planet around the sun. Note that units are in AU
    It takes in 7 parameters of those 6 are integers and 1 is a string. Returns 3 plots, 2 velocity plot
    and one plot describing the orbit(x vs y)"""
    # constants
    GM = 39.5
    mass = 1.65e-7 #the mass of the planet in solar
    x= [initx]
    y= [inity]
    vx= [initvx]
    vy= [initvy]
    radius = (initx**2+inity**2)**1.5
    radius2 = (initx**2+inity**2)
    rpow3=[radius]
    rpow2=[radius2]
    #for conservation of angular momentum
    initL = (initx*(mass*initvy)-inity*(mass*initvx))
    l=[initL]
    time=[0]
    delta_t = step  # the step time
    year = 1

    # plotting the graphs
    def plotting(x1, y1, vx1, vy1, l1, time1,titlexy,word,type):
        fig = plt.figure()
        ax = fig.add_subplot()
        if type == 0:
            fig1 = plt.figure(1)
            plt.plot(x1, y1)
            plt.plot(x1[len(times)],y1[len(times)],'ro',label = "perihelion") #plotting the furtherst point from the orbit
            plt.plot(x1[0],y1[0],'go',label="start of orbit")
            plt.plot(0,0,marker ='o',color ='gold',label ="the sun")
            plt.title(name + " orbiting the Sun for " + str(titlexy) +" "+ word)
            plt.xlabel("x (AU)")
            plt.ylabel("y (AU)")
            plt.axis('square')
            plt.legend(loc ="upper left", frameon ="False")
            plt.show()
            fig4 = plt.figure(4)
            plt.plot(time1, l1, '-b')
            plt.suptitle("Conservation of Angular momentum", y=0.99)
            plt.xlabel("Time")
            plt.ylabel("Angular momentum(kgAU/yr")
            plt.show()
            return fig1,fig4
        else:
            fig1 = plt.figure(1)
            plt.plot(x1, y1)
            #plt.plot(x1[len(times)], y1[len(times)], 'ro',label="perihelion")
            #plt.plot(x1[0], y1[0], 'go', label="start of orbit")
            plt.plot(0, 0, marker='o', color='gold', label="the sun")
            plt.title(name + " orbiting the Sun for " + str(titlexy) + " " + word)
            plt.xlabel("x (AU)")
            plt.ylabel("y (AU)")
            plt.axis('square')
            plt.legend(loc="upper left", frameon="False")
            plt.show()
            fig2 =plt.figure(2)
            plt.plot(time1, vx1)
            plt.title("X-component of Velocity vs Time")
            plt.xlabel("Time(yr)")
            plt.ylabel("Vx (AU/yr)")
            plt.show()
            fig3 =plt.figure(3)
            plt.plot(time1, vy1)
            plt.title("Y-component of Velocity vs Time")
            plt.xlabel("Time(yr)")
            plt.ylabel("Vy (AU/yr)")
            plt.show()
            # showing angular momentum is conserved
            fig4 =plt.figure(4)
            plt.plot(time1, l1,'-b')
            plt.suptitle("Conservation of Angular momentum", y=0.99)
            plt.xlabel("Time")
            plt.ylabel("Angular momentum(kgAU/yr")
            plt.show()
        return fig1,fig2,fig3,fig4
    if type == 0:
        year = 0
        alpha = 0.01
        for i in range(0,2):
            year+= 5
            times = np.arange(0,year,delta_t)
            print("GR time")
            #for different years we will use a loop
            for i in range(len(times)):
            #adds the new distance cubed
                rpow3.append((x[-1] ** 2 + y[-1] ** 2) ** 1.5)
            #adds the new distance sqaured
                rpow2.append((x[-1]**2+y[-1]**2))
            # adds the new x component of velocity
                vx.append((vx[-1] - GM *(1+alpha/rpow2[-1])*x[-1] * delta_t / rpow3[-1]))
            # adds the new x position of the planet
                x.append((x[-1] + vx[-1] * delta_t))
            # adds the new y component of velocity
                vy.append((vy[-1] - GM * (1+alpha/rpow2[-1])*y[-1] * delta_t / rpow3[-1]))
            # adds the new y position of the planet
                y.append((y[-1] + vy[-1] * delta_t))
            # adds the step to the previous time, results in current time.
                time.append(time[-1] + delta_t)
            # showing the angular momentum is conserved
                l.append((x[i+1] * mass * vy[i+1]) - (y[i+1] * mass * vx[i+1]))
            plotGR = plotting(x,y,vx,vy,l,time,year,"years",0)
    else:
        year = 1
        times = np.arange(0, year, delta_t)  # the number of iterations
        for i in range(len(times)):
            #adds new distance cubed
            rpow3.append((x[-1] ** 2 + y[-1] ** 2) ** 1.5)
            # adds the new x component of velocity
            vx.append((vx[-1]-GM*x[-1]*delta_t/rpow3[-1]))
            # adds the new x position of the planet
            x.append((x[-1]+vx[-1]*delta_t))
            # adds the new y component of velocity
            vy.append((vy[-1]-GM*y[-1]*delta_t/rpow3[-1]))
            # adds the new y position of the planet
            y.append((y[-1]+vy[-1]*delta_t))
            # adds the step to the previous time, results in current time.
            time.append(time[-1]+delta_t)
            #calculating the angular momentum
            l.append((x[-1]*mass*vy[-1])-(y[-1]*mass*vx[-1]))
        plotsR = plotting(x,y,vx,vy,l,time,year,"year",1)

#without GRcorrection part c
Mercury = orbit(0.47,0.0,0.0,8.17,0.0001,1,"Mercury")
#with GR correction part d
MercuryGR= orbit(0.47,0.0,0.0,8.17,0.0001,0,"Mercury")
