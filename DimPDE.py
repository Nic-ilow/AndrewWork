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
parser.add_argument('-p0','--p0',default=1)
parser.add_argument('-pscale','--pscale',default=.01)
parser.add_argument('-dx','--dx',default=7.0)
parser.add_argument('-D0','--D0',default = 2.7e5)
parser.add_argument('-xmax','--xmax',default=4.0e3)
parser.add_argument('-dt','--dt',default=1.0e-5)
parser.add_argument('-tmax','--tmax',default=32.0)
parser.add_argument('-arate','--arate',default=1.0)
parser.add_argument('-am','--am',default=13.0)
args = parser.parse_args()

p0 = float(args.p0)
dx = float(args.dx)
xmax = float(args.xmax)
dt = float(args.dt)
tmax = float(args.tmax)
D0 = float(args.D0)
arate = float(args.arate)
am = float(args.am)
pscale = float(args.pscale)


colorList = ['b','g','r','c','m','k','y']


x = np.arange(0,xmax,dx)
t = np.arange(0,tmax,dt)

width = np.size(x)

### pTot and aTot plotting array
tArray = np.logspace(-200,0,num=200,base=1.1)*tmax
### Initializing counters and a timer
tstart = time.time()
counter=0
### Array Pre-allocation
p = np.zeros(width+1)
p_1 = np.zeros_like(p)
       
telapsed = 0
print('p0hat =',p0, '\n dtbar = ',dt, '\n dxbar = ',dx, '\n pscale = ',pscale)
acetyl = np.zeros(width)
pTot = np.zeros(np.size(tArray))
aTot = np.zeros_like(pTot)
### Boundary Condition
p[0] = p0
p[width-1] = 0
counter2 = 0
SliceTimes = [2,30]
pArray = []
aArray = []
while telapsed<=tmax: # Iterating the system of tmax amount of seconds

        while (telapsed - tArray[counter2])>=0:
                pTot[counter2] = sum((p[1:width])) * dx
                aTot[counter2] = sum(acetyl[1:width]) * dx
                if counter2<np.size(tArray)-1:
                        counter2+=1
        #while any(telapsed-SliceTimes
        pscaler = 1+p/pscale+(p/pscale)**2
        p_1[1:width] = p[1:width] +   dt * D0 * ( ( (p[0:width-1] - p[1:width]) / ((pscaler[0:width-1]+pscaler[1:width])/2) + (p[2:width+1] - p[1:width]) / ((pscaler[2:width+1]+pscaler[1:width])/2 )) / (dx**2))
        acetyl[0:width] = acetyl[0:width] + (am - acetyl[0:width] ) * arate * dt * p[0:width] ### NO SFD

        telapsed += dt

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        p[width]=p[width-1] # Closed Tube Boundary
print((time.time()-tstart)/60)

counter = 0

plt.figure(1)
plt.plot(x/xmax,p[0:width])
plt.xlabel('x/L')
plt.ylabel('Density')

plt.figure(2)
plt.plot(x/xmax,acetyl/am)
plt.xlabel('x/L')
plt.ylabel('Acetylation')

plt.figure(3)
plt.loglog(tArray[:199] , pTot[:199])
plt.xlabel('t (s)')
plt.ylabel('Total Enzymes in Tubule')

plt.figure(4)
plt.loglog(tArray[:199] , aTot[:199])
plt.xlabel('t (s)')
plt.ylabel('Total Acetylated Sites ')

plt.show()
