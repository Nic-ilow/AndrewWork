import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis

parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-D0','--D0',default=2.7E5)
parser.add_argument('-dx','--dx',default=7.0)
parser.add_argument('-tmax','--tmax',default=2.0)
parser.add_argument('-width','--width',default=32)
parser.add_argument('-p0','--p0',default=1)
parser.add_argument('-dt','--dt',default=1E-6)
args = parser.parse_args()

### Relevant constants
width = int(args.width)
D0 = float(args.D0)
dx = float(args.dx)
tmax = float(args.tmax)
p0 = float(args.p0)
dt = float(args.dt)

D = D0*dt/(dx**2)
### Plotting Cutoffs
epsilon = dt/1000
t0 = 1E-1
t1 = 2E-1
t2 = 4E-1
t3 = 8E-1
t4 = 16E-1
### Concentration Variables
p = np.random.rand(width)
p_1 = np.zeros_like(p)
x=np.arange(width)
### Density Tracking
occupationDensity = np.zeros_like(p)

### Boundary Condition
p[0] = p0
p[width-1] = 0

telapsed = 0

def progress(p):

        newP1 = p[1] + ( D0 * ( p[0] - 2*p[1] + p[2] ) / (dx**2) * dt )
        return newP1

def rhopass(i):
        if i==width-1:
                pPass = np.array( [ p[i-1] , p[i] , p[i] ] )
        elif i==0:
                pPass = np.array( [ p[i] , p[i] , p[i+1] ] )
        else:
                pPass = np.array( [ p[i-1] , p[i] , p[i+1] ] )
        
        return pPass
np.vectorize(progress)

while telapsed<tmax:

        p_1[1:width-1] = p[1:width-1] + D * ( p[0:width-2] + p[2:width] - 2*p[1:width-1] ) 

        telapsed += dt
        if abs(telapsed-t0)<epsilon:
                plt.figure()
                plt.scatter(x,occupationDensity/telapsed,marker='o',label=('t= t0')) 
                plt.plot(x,scis.erfc(x/np.sqrt(4*D0*telapsed)))
                print(telapsed)
        if abs(telapsed-t1)<epsilon:
                plt.figure()
                plt.scatter(x,occupationDensity/telapsed,marker='o',label=('t= t1'))
                plt.plot(x,scis.erfc(x/np.sqrt(4*D0*telapsed))) 
                print(telapsed)
        if abs(telapsed-t2)<epsilon:
                plt.figure()
                plt.scatter(x,occupationDensity/telapsed,marker='o',label=('t= t2'))
                plt.plot(x,scis.erfc(x/np.sqrt(4*D0*telapsed)))
                print(telapsed)
        if abs(telapsed-t3)<epsilon:
                plt.figure()
                plt.scatter(x,occupationDensity/telapsed,marker='o',label=('t= t3'))
                plt.plot(x,scis.erfc(x/np.sqrt(4*D0*telapsed)))
                print(telapsed)
        if abs(telapsed-t4)<epsilon:
                plt.figure()
                plt.scatter(x,occupationDensity/telapsed,marker='o',label=('t= t4'))
                plt.plot(x,scis.erfc(x/np.sqrt(4*D0*telapsed)))
                print(telapsed)
        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        #p[width-1] = 0
        occupationDensity += p*dt

occupationDensity /= tmax
plt.figure()
plt.scatter(x,occupationDensity,marker='o',label='t= tmax')
plt.plot(x,scis.erfc(x/np.sqrt(4*D0*tmax)))
plt.legend()
plt.show()

