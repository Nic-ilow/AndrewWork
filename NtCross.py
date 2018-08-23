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
dxbarArray = np.logspace(-2,2,5,base=2)
p2 = np.logspace(-1,2,4,base=2) 
### Initializing counters and a timer 
tstart = time.time()
counter=0
### Array Pre-allocation
for dxbar in dxbarArray:
        counter+=1
        counter3 = 0
        
        dtbar = min(0.5/p0hat , 0.10*dxbar**2)
        print('p0hat = ',p0hat)
        print('dtbar = ',dtbar)
        print('dxbar = ',dxbar)
        epsilon = dtbar/2
        tArray = np.exp(np.arange(0,np.log(tbarmax),dtbar))-(1/dtbar-2)*dtbar
        xbar = np.arange(0,xbarmax,dxbar)
        tbar = np.arange(0,tbarmax,dtbar)
        width = np.size(xbar)       

        pTotTemp = np.zeros( (np.size(tArray) , np.size(p2)) )
        for p0hat in p2:
           
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
                                if counter2<np.size(tArray)-1:
                                        counter2+=1

                        pscaler = 1+p+p**2     
#                        p_1[1:width] = p[1:width] +   dtbar *  ( (p[0:width-1]/pscaler[0:width-1] + p[2:width+1]/pscaler[2:width+1] - 2*p[1:width]/pscaler[1:width]) / (dxbar**2))  
### NEW METHOD BELOW
                        p_1[1:width] = p[1:width] +   dtbar *  ( ( (p[0:width-1] - p[1:width]) / (pscaler[0:width-1]+pscaler[1:width])/2 + (p[2:width+1] - p[1:width]) / (pscaler[2:width+1]+pscaler[1:width])/2 ) / (dxbar**2)) 
                        acetyl[0:width] = acetyl[0:width] + (1 - acetyl[0:width] ) * dtbar * p[0:width] ### NO SFD
                         
                        telapsed += dtbar                     

                        #counter2+=1

                        p,p_1 = p_1,p # Updating concentration array
                        p[0] = p0hat # resetting the boundary condition
                        p[width]=p[width-1] # Closed Tube Boundary
 
                pTotTemp[:,counter3] = pTot
                '''
                if counter3==3:
                        idx = np.argwhere(np.diff(np.sign(pTotTemp[:,3] - pTotTemp[:,0]))).flatten()
                        plt.figure(10)
                        plt.scatter(dxbar,tArray[idx])
                        idx2 = np.argwhere(np.diff(np.sign(pTotTemp[:,3] - pTotTemp[:,1]))).flatten()
                        plt.figure(11)
                        plt.scatter(dxbar,tArray[idx2])
                '''                        
                plt.figure(counter)
                plt.scatter(xbar,p[0:width]/p0hat,s=2,label=('dxbar = %.2e p0hat = %.2e'%(dxbar,p0hat)))
                
                plt.figure(np.size(dxbarArray)+counter)
                plt.loglog(tArray , pTot,label=('dxbar = %.2e p0hat = %.2e'%(dxbar,p0hat)))
                              
                print((time.time()-tstart)/60)
                counter3 += 1

for i in range(np.size(dxbarArray)):
        plt.figure(i+1)
        plt.xlabel('Dimensionless length')
        plt.ylabel('phat/phat0')
        plt.legend()

        plt.figure(i+1+np.size(dxbarArray))
        plt.xlabel('Dimensionless Time')
        plt.ylabel('N(t)')
        plt.legend()

plt.show()

