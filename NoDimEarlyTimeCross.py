import math as m
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import scipy.integrate as scii
import scipy.optimize as scio
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
colorList = ['b','g','r','c','m','k','y']


xbar = np.arange(0,xbarmax,dxbar)
tbar = np.arange(0,tbarmax,dtbar)

width = np.size(xbar)       

### pTot and aTot plotting array
dtbar = min(0.5/p0hat , 0.10*dxbar**2)
tArray = np.exp(np.arange(0,np.log(tbarmax),dtbar))-(1/dtbar-2)*dtbar
dxbarArray = np.logspace(-1,-10,10,base=1.2)
p2 = np.array([.0625,64])
pArray = np.zeros((np.size(p2) ,np.size(dxbarArray), 2))        
### Initializing counters and a timer 
tstart = time.time()
counter=0
tcut = .05
for p0hat in p2:
        counter2 = 0
        for dxbar in dxbarArray:
                dtbar = min(0.5/p0hat , 0.10*dxbar**2)
                if dtbar==0.5/p0hat:
                        print('p0Hat Limit')
                else:
                        print('dxbar Limit')
                epsilon = dtbar/2
                p = np.zeros(width+1)
                p_1 = np.zeros_like(p)
                acetyl = np.zeros(width)
                
                telapsed = 0
                print('p0hat =',p0hat, '\n dtbar = ',dtbar, '\n dxbar = ',dxbar) 
                pTot = np.zeros(np.size(tArray))
                aTot = np.zeros_like(pTot)
                ### Boundary Condition
                p[0] = p0hat
                p[width-1] = 0
         
                while telapsed<=tcut: # Iterating the system of tmax amount of seconds

                        pscaler = 1+p+p**2     
                        p_1[1:width] = p[1:width] +   dtbar *  ( ( (p[0:width-1] - p[1:width]) / ((pscaler[0:width-1]+pscaler[1:width])/2) + (p[2:width+1] - p[1:width]) / ((pscaler[2:width+1]+pscaler[1:width])/2 )) / (dxbar**2))  
                        acetyl[0:width] = acetyl[0:width] + (1 - acetyl[0:width] ) * dtbar * p[0:width] ### NO SFD
                         
                        telapsed += dtbar                     

                        p,p_1 = p_1,p # Updating concentration array
                        p[0] = p0hat # resetting the boundary condition
                        p[width]=p[width-1] # Closed Tube Boundary
                temp = np.copy(p)
                pArray[counter,counter2,:] = temp[1:3]
                counter2+=1
        counter+=1

tend = time.time()
print((tend-tstart)/60)
for i, value1 in enumerate(p2):
        for j in range(2):
                plt.figure(1)
                plt.plot(dxbarArray,pArray[i,:,j],label=('position %d for p0hat = %.2f'%(j,value1)))
plt.figure(1)
plt.legend()
plt.xlabel('dxbar')
plt.ylabel('phat')
plt.title('t= %.2f'%tcut)
plt.show()
