import math as m
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import scipy.integrate as scii
import pickle
import time

parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-D','--D',default=2.7E-13)
parser.add_argument('-dx','--dx',default=7.0)
parser.add_argument('-tmax','--tmax',default=1.0)
parser.add_argument('-p0','--p0',default=1)
parser.add_argument('-arate','--arate',default=100)
parser.add_argument('-am','--am',default=1)
parser.add_argument('-pscale','--pscale',default=1000)
parser.add_argument('-dt','--dt',default=1e-6)
args = parser.parse_args()

global p0,arate,acetylmultiplicity,D
### Analytical Acetylation Function

def acetylAnalytical(x,telapsed):
        z = x / np.sqrt(4*D*telapsed)
        acetyl = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
        return  acetyl

def densityAnalytical(x,telapsed):
        z = x / np.sqrt(4*D*(telapsed))
        return scis.erfc(z)


### Relevant constants
D = float(args.D)
dx = float(args.dx)
tmax = float(args.tmax)
p0 = float(args.p0)
arate = float(args.arate)
acetylmultiplicity = int(args.am)
pscale = float(args.pscale)
am = float(args.am)

### Length and Step #### SCALING STUFF OUTs
L=3000e-9
dx = dx*1e-9
x=np.arange(0,L,dx)

dt = float(args.dt)
t=np.arange(0,tmax,dt)

tc = 1/(am*arate)
dtbar = dt/tc
tmaxbar = tmax/tc
tbar = t/tc

xc = np.sqrt(D*tc)  
xbar = x/xc
dxbar = dx/xc


#### SCALE PARAMETERS
p0hat = p0/(1+p0/pscale+(p0/pscale)**2)## Calculated values, D for phenomenological, and Stability limit for dt

epsilon = dt/2
width = np.size(xbar)       

### pTot and aTot plotting array
tArray = np.arange(1*dt,tmax+dt,dt*10)

### Initializing counters and a timer 
tstart = time.time()
counter=0
### Array Pre-allocation
p = np.zeros(width+1)
p_1 = np.zeros_like(p)
telapsed = 0

acetyl = np.zeros(width)
pTot = np.zeros((int(int(tmax/dt+1)/10+1)))
aTot = np.zeros_like(pTot)
aTotAnalytical =  np.zeros_like(aTot)
pTotAnalytical = np.zeros_like(aTotAnalytical)

### Boundary Condition
p[0] = p0hat
p[width-1] = 0
counter2 = 0

while telapsed<=tmaxbar: # Iterating the system of tmax amount of seconds
        pscaler = 1+p/pscale+(p/pscale)**2     
        p_1[1:width] = p[1:width] +   dtbar *  ( (p[0:width-1]*pscaler[0:width-1] + p[2:width+1]*pscaler[2:width+1] - 2*p[1:width]*pscaler[1:width]) / (dxbar**2))             ### NO SFD
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dtbar * p[0:width] ### NO SFD
         
        telapsed += dtbar                     

        #counter2+=1

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0hat # resetting the boundary condition
        p[width]=p[width-1] # Closed Tube Boundary

        if (counter2%10==0):
                pTot[int(counter2/10)] = sum(p[1:width]) * dxbar
                aTot[int(counter2/10)] = sum(acetyl) * dxbar
                
                z = xbar/np.sqrt(4*D*(telapsed))
                
                #pTotAnalytical[int(counter2/10)] = p0*np.sqrt(4*D*(telapsed)/np.pi)
                #pTotAnalytical[int(counter2/10)-1] = scii.quad(densityAnalytical,0,x[width-1],args=(telapsed),epsabs=1e-15)[0]
                #pTotAnalytical[int(counter2/10)-1] = sum(scis.erfc(z))
                #acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
                #aTotAnalytical[int(counter2/10)] = scii.quad(acetylAnalytical,0,xbar[width-1],args=(telapsed),epsabs=1e-15)[0] 
        counter2+=1
        #p_2[0] = p0
print((time.time()-tstart)/60)
#counter+=1
#print(telapsed)
#tend=time.time()
#trun = (tend-tstart)/60 #In minutes
#print(trun)

#print('residual')
#print(pTot2[0]-pTotAnalytical[0])
#print('max residual')
#print(max(abs(pTot2-pTotAnalytical)))
plt.scatter(xbar,p[0:width],s=1)
plt.show()
'''
### PLOTTING
parameterString = 'arate=%.2f \n L=%.2f \n dx=%.2f'%(arate,L,dx)
pTotResidual = np.zeros_like(dxArray)
#
plt.scatter(x/L , p[0:width] , marker=',', s=2 , label='SFD') 
plt.scatter(x/L , p_2[0:width] , marker='*' , s=2 , label='No SFD')
plt.xlabel('x/L')
plt.ylabel('Density')
plt.legend()

plt.figure()
plt.scatter(x/L , acetyl , marker=',' , s=2 , label='SFD')
plt.scatter(x/L , acetyl , marker='*' , s=2 , label='No SFD')
plt.xlabel('x/L')
plt.ylabel('Acetylation')
plt.legend()


plt.figure()
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.scatter(tArray , pTot , marker=',' , s=2 , label='SFD')
plt.scatter(tArray , pTot2 , marker='*' , s=2 , label='No SFD')
plt.xlabel('time (s)')
plt.ylabel('N(t)')
plt.legend()


plt.figure()
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.scatter(tArray , aTot , marker=',' , s=2 , label='SFD')
plt.scatter(tArray , aTot2 , marker='*' , s=2 , label='No SFD')
plt.xlabel('time (s)')
plt.ylabel('A(t)')
plt.legend()

#pTotRatio = 
plt.figure()
ax=plt.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
plt.scatter(tArray , pTot/pTot2 , marker=',' , s=2 , label='N(t)')
plt.scatter(tArray , aTot/aTot2 , marker='*' , s=2 , label='A(t)')
plt.xlabel('time (s)')
plt.ylabel('Ratio')
plt.title('Ratio of SFD vs No SFD')
plt.legend()

plt.figure()
ax=plt.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
#plt.scatter(tArray , pTot/pTotAnalytical , marker=',' , s=2 , label='SFD over Analytical')
plt.scatter(tArray , abs(pTot2-pTotAnalytical) , marker='*' , s=2 , label='No SFD over Analytical')
plt.xlabel('time (s)')
plt.ylabel('Ratio')
plt.title('Ratio of SFD and No SFD vs Analytical')
plt.legend()


plt.show()
'''
