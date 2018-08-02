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
parser.add_argument('-xbarmax','--xbarmax',default=1800.0)
parser.add_argument('-dtbar','--dtbar',default=.1)
parser.add_argument('-tbarmax','--tbarmax',default=1.0e5)
args = parser.parse_args()

p0hat = float(args.p0hat)
dxbar = float(args.dxbar)
xbarmax = float(args.xbarmax)
dtbar = float(args.dtbar)
tbarmax = float(args.tbarmax)
epsilon = dtbar/2

### pTot and aTot plotting array
tArray = np.exp(np.arange(0,np.log(tbarmax),dtbar))
dxbarArray = np.logspace(-5,-1,5,base=2)

### Initializing counters and a timer 
tstart = time.time()
counter=0
### Array Pre-allocation
plt.figure(3)
ax3 = plt.gca()
plt.figure(4)
ax4 = plt.gca()

for dxbar in dxbarArray: 

        xbar = np.arange(0,xbarmax,dxbar)
        tbar = np.arange(0,tbarmax,dtbar)
        width = np.size(xbar)       
 
        p = np.zeros(width+1)
        p_1 = np.zeros_like(p)
        telapsed = 0

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
                        #z = xbar/np.sqrt(4*D*(telapsed))
                        
                        #pTotAnalytical[int(counter2/10)] = p0*np.sqrt(4*D*(telapsed)/np.pi)
                        #pTotAnalytical[int(counter2/10)-1] = scii.quad(densityAnalytical,0,x[width-1],args=(telapsed),epsabs=1e-15)[0]
                        #pTotAnalytical[int(counter2/10)-1] = sum(scis.erfc(z))
                        #acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
                        #aTotAnalytical[int(counter2/10)] = scii.quad(acetylAnalytical,0,xbar[width-1],args=(telapsed),epsabs=1e-15)[0] 
                        if counter2<np.size(tArray)-1:
                                counter2+=1

                pscaler = 1+p+p**2     
                p_1[1:width] = p[1:width] +   dtbar *  ( (p[0:width-1]/pscaler[0:width-1] + p[2:width+1]/pscaler[2:width+1] - 2*p[1:width]/pscaler[1:width]) / (dxbar**2))  
                acetyl[0:width] = acetyl[0:width] + (1 - acetyl[0:width] ) * dtbar * p[0:width] ### NO SFD
                 
                telapsed += dtbar                     

                #counter2+=1

                p,p_1 = p_1,p # Updating concentration array
                p[0] = p0hat # resetting the boundary condition
                p[width]=p[width-1] # Closed Tube Boundary

        plt.figure(1)
        plt.scatter(xbar,p[0:width]/p0hat,s=1,label=('dxbar = %.2e'%dxbar))
        plt.xlabel('Dimensionless length')
        plt.ylabel('phat/phat0')

        plt.figure(2)
        plt.scatter(xbar,acetyl,s=1,label=('dxbar = %.2e'%dxbar))
        plt.xlabel('Dimensionless length')
        plt.ylabel('ahat')
        
        plt.figure(3)
        plt.scatter(tArray , pTot,s=1,label=('dxbar = %.2e'%dxbar))
        plt.xlabel('Dimensionless Time')
        plt.ylabel('N(t)')
       
        plt.figure(4)
        plt.scatter(tArray , aTot,s=1,label=('dxbar = %.2e'%dxbar))
        plt.xlabel('Dimensionless Time')
        plt.ylabel('A(t)')
        
        counter+=1
        print((time.time()-tstart)/60)

plt.figure(1)
plt.legend()
plt.figure(2)
plt.legend()
plt.figure(3)
ax3.set_xscale('log')
ax3.set_yscale('log')
plt.legend()
plt.figure(4)
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.legend()

plt.show()
