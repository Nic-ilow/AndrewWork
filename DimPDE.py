import math as m
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import scipy.integrate as scii
import scipy.optimize as scio
import pickle
import datetime
import time
import os
import shutil

now = datetime.datetime.now()
dateFormat = 'Sim_{0}_{1}_{2}_{3}_{4}'.format(now.day,now.month,now.year,now.hour,now.minute)

if os.path.exists(dateFormat):
        shutil.rmtree(dateFormat)
os.makedirs(dateFormat)

parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-p0','--p0',default=1)
#parser.add_argument('-pscale','--pscale',default=.01)
parser.add_argument('-kon','--kon',default=1000)
parser.add_argument('-koff','--koff',default=10)
parser.add_argument('-dx','--dx',default=7.0)
parser.add_argument('-D0','--D0',default = 2.7e5)
parser.add_argument('-xmax','--xmax',default=3584.0)
parser.add_argument('-dt','--dt',default=1.0e-5)
parser.add_argument('-tmax','--tmax',default=1.0)
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
#pscale = float(args.pscale)
kon = float(args.kon)
koff = float(args.koff)
if kon==0:
        pscale = 10000000
else:
        Kd = koff/kon
        fractionfree = Kd/(1+Kd)

        khop = 2.0*D0/(dx*dx)/fractionfree
        pscale = np.sqrt((kon+koff)/khop)*(1+Kd)
print('The effective pscale is:',pscale)
colorList = ['b','g','r','c','m','k','y']


x = np.arange(0,xmax,dx)
t = np.arange(0,tmax,dt)

width = np.size(x)

simInfo = [dx,width,kon,koff,dt,arate,tmax]

### pTot and aTot plotting array
plotCuts = [.1, .5, 1, 2, 4, 8, 16, 30]
netCuts = list(np.logspace(-200,0,num=200,base=1.1)*tmax)
### Initializing counters and a timer
tstart = time.time()
counter=0
### Array Pre-allocation
p = np.zeros(width+1)
p_1 = np.zeros_like(p)
       
telapsed = 0
print('p0hat =',p0, '\n dtbar = ',dt, '\n dxbar = ',dx, '\n pscale = ',pscale)
acetyl = np.zeros(width)
NTot = []
ATot = []
### Boundary Condition
p[0] = p0
p[width-1] = 0
counter2 = 0
pArray = []
aArray = []
plotPoints = []
totPoints = []
while telapsed<=tmax: # Iterating the system of tmax amount of seconds
        if len(plotCuts) > 0:
                while (telapsed + dt - plotCuts[0]) > 0:
                        pArray.append(np.copy(p[0:width]))
                        aArray.append(np.copy(acetyl))
                        plotPoints.append(plotCuts.pop(0))
                        if len(plotCuts) == 0:
                                break
        if len(netCuts) > 0:
                while (telapsed+dt-netCuts[0])>0:
                        NTot.append(sum(np.copy(p[1:width])))
                        ATot.append(sum(np.copy(acetyl[1:])))
                        totPoints.append(netCuts.pop(0))
                        if len(netCuts)==0:
                                break
        #while any(telapsed-SliceTimes
        pscaler = 1+p/pscale+(p/pscale)**2
        p_1[1:width] = p[1:width] +   dt * D0 * ( ( (p[0:width-1] - p[1:width]) / ((pscaler[0:width-1]+pscaler[1:width])/2) + (p[2:width+1] - p[1:width]) / ((pscaler[2:width+1]+pscaler[1:width])/2 )) / (dx**2))
        acetyl[0:width] = acetyl[0:width] + (am - acetyl[0:width] ) * arate * dt * p[0:width] ### NO SFD

        telapsed += dt

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        p[width]=p[width-1] # Closed Tube Boundary
print((time.time()-tstart)/60)

np.save(os.path.join(dateFormat, 'DENS'), pArray)
np.save(os.path.join(dateFormat, 'ACETYL'), aArray)
np.save(os.path.join(dateFormat, 'NTOT'), NTot)
np.save(os.path.join(dateFormat, 'ATOT'), ATot)
np.savetxt(os.path.join(dateFormat,'SIMINFO'),simInfo,delimiter=',')
if os.path.exists('ContinuousData'):
        shutil.move(dateFormat,'ContinuousData')
'''
def aAnalytic(x,Beta,t):
        global p0,arate,D0
        z = x/np.sqrt(4*D0*t)
        afit = 1 - np.exp(-p0*Beta*arate*t*( ((1+2*z*z)*scis.erfc(z)) - ((2*z*np.exp(-(z*z)))/np.sqrt(np.pi)) ) )
        return afit

def pAnalytic(x,t):
        z = x/np.sqrt(4*D0*t)
        pfit = scis.erfc(z) * p0
        return pfit

pfit = pAnalytic(x,tmax)
afit = aAnalytic(x,1,tmax)
plt.figure(1)
plt.scatter(x/xmax,p[0:width],s=2,label=('Finite Difference Model at t=%.2f s')%(tmax))
plt.plot(x/xmax,pfit,c='r',ls='dashed',lw=4,label='Simple Diffusion Fit')
plt.xlabel('x/L')
plt.ylabel('Density')
plt.legend()

plt.figure(2)
plt.scatter(x/xmax,acetyl/am,s=2,label=('Finite Difference Model at t=%.2f s'%(tmax)))
plt.plot(x/xmax,afit,c='r',ls='dashed',lw=4,label='Simple Diffusion Fit')
plt.xlabel('x/L')
plt.ylabel('Acetylation')
plt.legend()
NtotAnalytic = p0*np.sqrt(4*D0/(dx*dx)*tArray/np.pi)
plt.figure(3)
plt.loglog(tArray[:199] , pTot[:199],label=('Finite Difference Model at t=%.2f s'%(tmax)))
plt.loglog(tArray[:199],NtotAnalytic[:199],label='Simple Diffusion Fit')
plt.xlabel('t (s)')
plt.ylabel('Total Enzymes in Tubule')
plt.legend()
ATotFit = []
for t in tArray:
        ATotFit.append(sum(aAnalytic(x,1,t)))
plt.figure(4)
plt.loglog(tArray[:199] , aTot[:199],label=('Finite Difference Model at t=%.2f s'%(tmax)))
plt.loglog(tArray[:199],ATotFit[:199],label='Simple Diffusion Fit')
plt.xlabel('t (s)')
plt.ylabel('Total Acetylated Sites ')
plt.legend()

plt.show()
'''
