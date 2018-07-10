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
L = 3000

### Calculated values, D for phenomenological, and Stability limit for dt
D = D0/(dx**2)
dt = (dx**2)/D0*.5
epsilon = dt/2
x = np.arange(0,L,dx)
width = np.size(x)       

### Plotting Cutoffs
base = 2
Samples = 6
tmaxArray = np.logspace(-1,Samples,Samples+2,base=base)

### Array Pre-allocation
p = np.zeros(width+1)
p_1 = np.zeros_like(p)
acetyl = np.zeros(width)
pTot = np.zeros(int(int(tmax/dt+1)/10)+1)
aTot = np.zeros_like(pTot)
pTotAnalytical = np.zeros_like(pTot)
aTotAnalytical = np.zeros_like(aTot)
densityExperiment = np.zeros( (np.size(tmaxArray) , width) )
acetylExperiment = np.zeros_like(densityExperiment)
densityResidual = np.zeros_like(densityExperiment)
acetylResidual = np.zeros_like(acetylExperiment)
### Boundary Condition
p[0] = p0
p[width-1] = 0

### Initializing 
telapsed = 0
print(telapsed)

### pTot and aTot plotting array
tArray = np.arange(0,tmax,dt*10)

### Initializing counters and a timer 
counter = 0
counter2 = 0
tstart = time.time()


while telapsed<=tmax: # Iterating the system of tmax amount of seconds

        Dsfd = (D/(1+(p/pscale)+((p/pscale)**2))) # Updating the Single File Diffusion effects based on density at positions
        
        p_1[1:width] = p[1:width] +  dt * ( Dsfd[0:width-1]*p[0:width-1] + Dsfd[2:width+1]*p[2:width+1] - Dsfd[1:width]*2*p[1:width]) 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

        telapsed += dt                     

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        p[width]=p[width-1] # Closed Tube Boundary
       
        if (counter%10)==0:
                pTot[int(counter/10)] = sum(p[0:width])
                pTotAnalytical[int(counter/10)] = p0 * np.sqrt(4*D*telapsed/np.pi)
        
                aTot[int(counter/10)] = sum(acetyl)
                aTotAnalytical[int(counter/10)] = scii.quad(acetylAnalytical,0,x[width-1],args=(telapsed))[0] 
        counter+=1

        if counter2<7:
                if abs(telapsed-tmaxArray[counter2])<=epsilon:
                        
                        plt.figure(1)
                        plt.plot(x/L,p[0:width],label='t=%.2f (s)'%tmaxArray[counter2])
                        
                        plt.figure(2)
                        plt.plot(x/L,acetyl,label='t=%.2f (s)'%tmaxArray[counter2])
                        
                        densityExperiment[counter2,:] = p[0:width]
                        acetylExperiment[counter2,:] = acetyl

                        z = x/np.sqrt(4*D0*telapsed) #Variable to make next equations easier to write                            
                        acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ) 
                        densityfit = scis.erfc(z)
                        
                        densityResidual[counter2,:] = densityExperiment[counter2,:] - densityfit
                        acetylResidual[counter2,:] = acetylExperiment[counter2,:] - acetylationfit
                        
                        print(telapsed)
                        counter2+=1

tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)

### PLOTTING
parameterString = ' dt (s)=%.2e \np0=%.2f \nLength of Tubule (nm)=%.2f  '%(dt,p0,L)


### DENSITY AT GIVEN TIMES
plt.figure(1)
plt.xlabel('x scaled (x/L)')
plt.ylabel('Density')
plt.legend()


### ACETYLATION AT GIVEN TIMES
plt.figure(2)
plt.xlabel('x scaled (x/L)')
plt.ylabel('Acetylation')
plt.legend()


### PLOTTING TOTAL DENSITY
plt.figure()
plt.xlabel('time (s)')
plt.ylabel('N(t)')
plt.plot(tArray,pTot)
plt.plot(tArray,pTotAnalytical)

### PLOTTING TOTAL ACETYLATION
plt.figure()
plt.xlabel('time (s)')
plt.ylabel('A(t)')
plt.plot(tArray,aTot)
plt.plot(tArray,aTotAnalytical)

plt.show()

'''
plt.figure()
ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.scatter(x/L,DensityResidual,s=1,label=('Residual at t=%.2f'%tmax))
plt.xlabel('xscaled (x/L)')
plt.ylabel('Residual of Density')
plt.legend()


plt.figure()
ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.scatter(x/L,AcetylResidual,s=1,label=('Residual at t=%.2f'%tmax))
plt.xlabel('x scaled (x/L)')
plt.ylabel('Residual of Acetylation')
plt.legend()

fig = plt.figure()
ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.scatter(tArray,NetDensityResidual,s=1)
plt.xlabel('time (s)')
plt.ylabel('Residual of Total Density')

plt.figure()
ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.scatter(tArray,NetAcetylResidual,s=1)
plt.xlabel('time (s)')
plt.ylabel('Max Residual of Total Acetylation')


plt.show()
'''
