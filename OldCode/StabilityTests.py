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
parser.add_argument('-dtbar','--dtbar',default=10)
parser.add_argument('-tbarmax','--tbarmax',default=1.0e5)
parser.add_argument('-a','--a',default=-4)
parser.add_argument('-b','--b',default=-4)
args = parser.parse_args()

a = int(args.a)
b = int(args.b)
p0hat = float(args.p0hat)
dxbar = float(args.dxbar)
xbarmax = float(args.xbarmax)
dtbar = float(args.dtbar)
tbarmax = float(args.tbarmax)

xbar = np.arange(0,xbarmax,dxbar)
tbar = np.arange(0,tbarmax,dtbar)

epsilon = dtbar/2
width = np.size(xbar)       

### pTot and aTot plotting array
tArray = np.exp(np.arange(0,np.log(tbarmax),dtbar))
p2 = np.logspace(a,b,b-a+1,base=2)
pArray = np.zeros((np.size(p2) , width))        
aArray = np.zeros_like(pArray)
pTotArray = np.zeros( ( np.size(p2) , np.size(tArray) ) )
aTotArray = np.zeros_like(pTotArray)
### Initializing counters and a timer 
tstart = time.time()
counter=0
### Array Pre-allocation
for p0hat in p2:
        p = np.zeros(width+1)
        p_1 = np.zeros_like(p)
        telapsed = 0
        dtbar = min(0.5/p0hat , 0.10*dxbar**2) 
        acetyl = np.zeros(width)
        pTot = np.zeros(np.size(tArray))
        aTot = np.zeros_like(pTot)
        ### Boundary Condition
        p[0] = p0hat
        p[width-1] = 0
        counter2 = 0
        print(dtbar)
       
        while telapsed<=tbarmax: # Iterating the system of tmax amount of seconds
                
                pscaler = 1+p+p**2     
                p_1[1:width] = p[1:width] +   dtbar *  ( (p[0:width-1]/pscaler[0:width-1] + p[2:width+1]/pscaler[2:width+1] - 2*p[1:width]/pscaler[1:width]) / (dxbar**2))  
                acetyl[0:width] = acetyl[0:width] + (1 - acetyl[0:width] ) * dtbar * p[0:width] ### NO SFD
                 
                telapsed += dtbar                     

                p,p_1 = p_1,p # Updating concentration array
                p[0] = p0hat # resetting the boundary condition
                p[width]=p[width-1] # Closed Tube Boundary
       
        pArray[counter,:] = p[0:width]
        aArray[counter,:] = acetyl[0:width]

        counter+=1
        print((time.time()-tstart)/60)

for i , value in enumerate(p2): 
        plt.figure(1)
        plt.scatter(xbar,pArray[i,0:width]/value,s=1,label=('p0hat = %.2e'%value))
        plt.xlabel('Dimensionless length')
        plt.ylabel('phat/phat0')

        plt.figure(2)
        plt.scatter(xbar,aArray[i,:],s=1,label=('p0hat = %.2e'%value))
        plt.xlabel('Dimensionless length')
        plt.ylabel('ahat')
      

plt.figure(1)
plt.legend()
plt.figure(2)
plt.legend()
plt.show()
