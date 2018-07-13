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


### pTot and aTot plotting array
tArray = np.arange(dt,tmax+dt,dt*10)

pscaleArray=np.logspace(-2,3,6,base=10)
#pscaleArray = np.linspace(0.1,1,9)
pTotArray = np.zeros( (np.size(pscaleArray) ,int(int(tmax/dt+1)/10+1) )) 
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
        pTot = np.zeros((int(int(tmax/dt+1)/10+1)))
        aTot = np.zeros_like(pTot)
        aTotAnalytical =  np.zeros_like(aTot)
        pTotAnalytical = np.zeros_like(aTotAnalytical)
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
               
                if (counter2%10==0):
                        pTot[int(counter2/10)] = sum(p[0:width]) 
                        aTot[int(counter2/10)] = sum(acetyl)
                        
                        if pscale==pscaleArray[-1]:
                                z = x/np.sqrt(4*D0*telapsed)
                                #pTotAnalytical[int(counter2/10)] = sum(scis.erfc(z))
                                pTotAnalytical[int(counter2/10)] = p0*np.sqrt(4*D0*telapsed/np.pi)/dx 
                                acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
                                #aTotAnalytical[int(counter2/10)] = sum(acetylationfit)
                                aTotAnalytical[int(counter2/10)] = scii.quad(acetylAnalytical,0,x[width-1],args=(telapsed),epsabs=1e-15)[0] 
                counter2+=1

        pTotArray[counter,:] = pTot
        aTotArray[counter,:] = aTot
        densityArray[counter,:] = p[0:width]
        acetylArray[counter,:] = acetyl
        
        print(telapsed)
        counter+=1

z = x/np.sqrt(4*D0*telapsed) #Variable to make next equations easier to write                                    

acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ) 
densityfit = scis.erfc(z)

#pTotAnalytical = p0 * np.sqrt((4*D*(tArray)/np.pi))

tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)

### PLOTTING
parameterString = 'arate=%.2f \n L=%.2f \n dx=%.2f'%(arate,L,dx)
for i in np.arange(np.size(pscaleArray)):
                
        plt.figure(1)
        plt.xlabel('x/L')
        plt.ylabel('p(x)')
        plt.scatter(x/L,densityArray[i,:],s=1,label=('pscale = %.2e'%pscaleArray[i]))

        plt.figure(2)
        plt.xlabel('x/L')
        plt.ylabel('a(x)')
        plt.scatter(x/L,acetylArray[i,:],s=1,label=('pscale = %.2e'%pscaleArray[i]))
        #plt.plot(x/L,acetylationfit)

        plt.figure(3)
        plt.xlabel('time (s)')
        plt.ylabel('N(t)')
        plt.scatter(tArray,pTotArray[i,:],s=1,label=('pscale = %.2e'%pscaleArray[i]))
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')


        plt.figure(4)
        plt.xlabel('time (s)')
        plt.ylabel('A(t)')
        plt.scatter(tArray,aTotArray[i,:],s=1,label=('pscale = %.2e'%pscaleArray[i]))
        ax=plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')

plt.figure(1)
plt.plot(x/L,densityfit,label='No SFD Analytical')
plt.legend()
plt.title(parameterString)

plt.figure(2)
plt.plot(x/L,acetylationfit,label='No SFD Analytical')
plt.legend()

plt.figure(3)
plt.plot(tArray,pTotAnalytical,label='No SFD Analytical')
plt.legend()

plt.figure(4)
plt.plot(tArray,aTotAnalytical,label='No SFD Analytical')
plt.legend()

pTotResidual = abs(pTotArray[i,:] / pTotAnalytical)
aTotResidual = abs(aTotArray[i,:] / aTotAnalytical)

plt.figure(5)
ax1=plt.gca()
plt.plot(tArray,(pTotResidual),marker='.',ms=3,label='Nt/NtAnalytical')
plt.plot(tArray,(aTotResidual),marker='*',ms=3,label='At/AtAnalytical')
ax1.set_xscale('log')
#ax1.set_yscale('log')
plt.xlabel('time (s)')
plt.ylabel('Ratio')
plt.legend()

plt.show()
