import math as m
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import scipy.integrate as scii
import pickle
import time

parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-D0','--D0',default=2.7E5)
parser.add_argument('-dx','--dx',default=7.0)
parser.add_argument('-tmax','--tmax',default=1.0)
parser.add_argument('-p0','--p0',default=1)
parser.add_argument('-arate','--arate',default=100)
parser.add_argument('-am','--am',default=1)
parser.add_argument('-pscale','--pscale',default=1000)
args = parser.parse_args()

global p0,arate,acetylmultiplicity,D0

### Relevant constants
D0 = float(args.D0)
dx = float(args.dx)
tmax = float(args.tmax)
p0 = float(args.p0)
arate = float(args.arate)
acetylmultiplicity = int(args.am)
pscale = float(args.pscale)

### Analytical Acetylation Function

def acetylAnalytical(x,telapsed):
        z = (x*dx) / np.sqrt(4*D0*telapsed)
        acetyl = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
        return  acetyl
       

### Length and Steps
L = 6000

### Calculated values, D for phenomenological, and Stability limit for dt
D = D0/(dx**2)
dt = (dx**2)/D0*.5
epsilon = dt/2
x = np.arange(0,L,dx)
width = np.size(x)       


### pTot and aTot plotting array
tArray = np.arange(0,tmax,dt*10)

pscaleArray=np.logspace(-5,2,8,base=2)
#pscaleArray = np.linspace(0.1,1,9)
pTotArray = np.zeros( (np.size(pscaleArray) , int(int(tmax/dt+1)/10)+1) ) 
aTotArray = np.zeros_like(pTotArray)
densityArray = np.zeros( (np.size(pscaleArray) , width) )
acetylArray = np.zeros_like(densityArray)

### Initializing counters and a timer 
counter = 0
tstart = time.time()

for pscale in pscaleArray:
        
        telapsed = 0
        ### Array Pre-allocation
        p = np.zeros(width+1)
        p_1 = np.zeros_like(p)
        acetyl = np.zeros(width)
        pTot = np.zeros(int(int(tmax/dt+1)/10)+1)
        aTot = np.zeros_like(pTot)
        aTotAnalytical =  np.zeros_like(aTot)
        ### Boundary Condition
        p[0] = p0
        p[width-1] = 0

        counter2 = 0
        while telapsed<=tmax: # Iterating the system of tmax amount of seconds

                Dsfd = (D/(1+(p/pscale)+((p/pscale)**2))) # Updating the Single File Diffusion effects based on density at positions
                
                p_1[1:width] = p[1:width] +  dt * ( Dsfd[0:width-1]*p[0:width-1] + Dsfd[2:width+1]*p[2:width+1] - Dsfd[1:width]*2*p[1:width]) 
                acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

                telapsed += dt                     

                p,p_1 = p_1,p # Updating concentration array
                p[0] = p0 # resetting the boundary condition
                p[width]=p[width-1] # Closed Tube Boundary
               
                if (counter2%10)==0:
                        pTot[int(counter2/10)] = sum(p[0:width]) 
                        aTot[int(counter2/10)] = sum(acetyl)
                        
                        if pscale==pscaleArray[-1]:
                                aTotAnalytical[int(counter2/10)] = scii.quad(acetylAnalytical,0,x[width-1],args=(telapsed+dt))[0] 
                counter2+=1

        pTotArray[counter,:] = pTot
        aTotArray[counter,:] = aTot
        densityArray[counter,:] = p[0:width]
        acetylArray[counter,:] = acetyl
        
        print(telapsed)
        counter+=1

z = x/np.sqrt(4*D0*telapsed+dt) #Variable to make next equations easier to write                                    

acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ) 
densityfit = scis.erfc(z)

pTotAnalytical = p0 * np.sqrt(4*D*(tArray+dt)/np.pi)

tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)

### PLOTTING
parameterString = ' dt (s)=%.2e \np0=%.2f \nLength of Tubule (nm)=%.2f  '%(dt,p0,L)

plt.figure(1)
plt.xlabel('x/L')
plt.ylabel('p(x)')
plt.plot(x/L,densityArray[0,:])
plt.plot(x/L,densityArray[1,:])
plt.plot(x/L,densityArray[2,:])
plt.plot(x/L,densityArray[3,:])
plt.plot(x/L,densityArray[4,:])
#plt.plot(x/L,densityfit)

plt.figure(2)
plt.xlabel('x/L')
plt.ylabel('a(x)')
plt.plot(x/L,acetylArray[0,:])
plt.plot(x/L,acetylArray[1,:])
plt.plot(x/L,acetylArray[2,:])
plt.plot(x/L,acetylArray[3,:])
plt.plot(x/L,acetylArray[4,:])
#plt.plot(x/L,acetylationfit)

plt.figure(3)
plt.xlabel('time (s)')
plt.ylabel('N(t)')
plt.plot(tArray,pTotArray[0,:])
plt.plot(tArray,pTotArray[1,:])
plt.plot(tArray,pTotArray[2,:])
plt.plot(tArray,pTotArray[3,:])
plt.plot(tArray,pTotArray[4,:])
#plt.plot(tArray,pTotAnalytical)
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')


plt.figure(4)
plt.xlabel('time (s)')
plt.ylabel('A(t)')
plt.plot(tArray,aTotArray[0,:])
plt.plot(tArray,aTotArray[1,:])
plt.plot(tArray,aTotArray[2,:])
plt.plot(tArray,aTotArray[3,:])
plt.plot(tArray,aTotArray[4,:])
#plt.plot(tArray,aTotAnalytical)
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()
