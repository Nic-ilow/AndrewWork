import math as m
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import scipy.integrate as scii
import pickle
import time

parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-p0hat','--p0hat',default=1e-4)
parser.add_argument('-dxbar','--dxbar',default=6.0)
parser.add_argument('-xbarmax','--xbarmax',default=1250.0)
parser.add_argument('-dtbar','--dtbar',default=.1)
parser.add_argument('-tbarmax','--tbarmax',default=1.0e4)
args = parser.parse_args()

p0hat = float(args.p0hat)
dxbar = float(args.dxbar)
xbarmax = float(args.xbarmax)
dtbar = float(args.dtbar)
tbarmax = float(args.tbarmax)

### Initializing counters and a timer 
tstart = time.time()
counter=0
### Array Pre-allocation
dxbarArray = np.logspace(1,-3,5,base=2)
for dxbar in dxbarArray:
        dtbar = min(0.5/p0hat , 0.10*dxbar**2)
        epsilon = dtbar/2
        tArray = np.exp(np.arange(0,np.log(tbarmax),dtbar))
        xbar = np.arange(0,xbarmax,dxbar)
        tbar = np.arange(0,tbarmax,dtbar)
        width = np.size(xbar)       
        p = np.zeros(width+1)
        p_1 = np.zeros_like(p)
        telapsed = 0
        print('p0hat =',p0hat, '\n dtbar = ',dtbar, '\n dxbar = ',dxbar) 
        acetyl = np.zeros(width)
        pTot = np.zeros(np.size(tArray))
        aTot = np.zeros_like(pTot)
        ### Boundary Condition
        p[0] = p0hat
        p[width-1] = 0
        counter2 = 0
       
        while telapsed<=tbarmax: # Iterating the system of tmax amount of seconds

                if abs(tArray[counter2]-telapsed)<=epsilon:
                        pTot[counter2] = sum((p[0:width])) * dxbar
                        aTot[counter2] = sum(acetyl) * dxbar

                        if counter2<np.size(tArray)-1:
                                counter2+=1

                pscaler = 1+p+p**2     
                p_1[1:width] = p[1:width] +   dtbar *  ( ( (p[0:width-1] - p[1:width]) / ((pscaler[0:width-1]+pscaler[1:width])/2) + (p[2:width+1] - p[1:width]) / ((pscaler[2:width+1]+pscaler[1:width])/2) ) / (dxbar**2))  
                acetyl[0:width] = acetyl[0:width] + (1 - acetyl[0:width] ) * dtbar * p[0:width] ### NO SFD
                 
                telapsed += dtbar                     

                p,p_1 = p_1,p # Updating concentration array
                p[0] = p0hat # resetting the boundary condition
                p[width]=p[width-1] # Closed Tube Boundary
        
        plt.figure(1)
        plt.scatter(xbar,p[0:width]/p0hat,s=2,label=('dxbar = %.2e'%dxbar))
        
        plt.figure(2)
        plt.scatter(xbar,acetyl,s=2,label=('dxbar = %.2e'%dxbar))
        
        plt.figure(3)
        plt.scatter(tArray , pTot,s=2,label=('dxbar = %.2e'%dxbar))

        plt.figure(4)
        plt.scatter(tArray , aTot,s=2,label=('dxbar = %.2e'%dxbar))
 
        counter+=1
        
        print((time.time()-tstart)/60)
       
plt.figure(1)
plt.legend()
plt.xlabel('Dimensionless length')
plt.ylabel('phat/phat0')
plt.title('p0hat = %.1e'%p0hat)


plt.figure(2)
plt.legend()
plt.xlabel('Dimensionless length')
plt.ylabel('ahat')
plt.title('p0hat = %.1e'%p0hat)

plt.figure(3)
ax3 = plt.gca()
ax3.set_xscale('log')
ax3.set_yscale('log')
plt.legend()
plt.xlabel('Dimensionless Time')
plt.ylabel('N(t)')
plt.title('p0hat = %.1e'%p0hat)

plt.figure(4)
ax4 = plt.gca()
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.legend()
plt.xlabel('Dimensionless Time')
plt.ylabel('A(t)')
plt.title('p0hat = %.1e'%p0hat)

plt.show()
