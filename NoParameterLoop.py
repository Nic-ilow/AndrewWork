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
parser.add_argument('-dt','--dt',default=1e-6)
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
        z = x / np.sqrt(4*D0*telapsed)
        acetyl = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
        return  acetyl

def densityAnalytical(x,telapsed):
        z = x / np.sqrt(4*D0*(telapsed))
        return scis.erfc(z)

### Length and Steps
L = 3000

### Calculated values, D for phenomenological, and Stability limit for dt
#dt = 1e-6
dt = float(args.dt)
epsilon = dt/2
print(D0*dt)
### Setting Tube Array Discretization
x = np.arange(0,L,dx)
width = np.size(x)       

### pTot and aTot plotting array
tArray = np.arange(1*dt,tmax+dt,dt*10)

### Initializing counters and a timer 
tstart = time.time()

telapsed = 0

### Array Pre-allocation
p = np.zeros(width+1)
p_1 = np.zeros_like(p)
p_2 = np.zeros(width+1)
p_3 = np.zeros_like(p_2)


acetyl = np.zeros(width)
acetyl2 = np.zeros(width)
#pTot = np.zeros(2)
pTot = np.zeros((int(int(tmax/dt+1)/10+1)))
aTot = np.zeros_like(pTot)
pTot2 = np.zeros_like(pTot)
aTot2 = np.zeros_like(pTot)
aTotAnalytical =  np.zeros_like(aTot)
pTotAnalytical = np.zeros_like(aTotAnalytical)


### Boundary Condition
p[0] = p0
p_2[0] = p0
p[width-1] = 0
counter2 = 0
while telapsed<=tmax: # Iterating the system of tmax amount of seconds

        Dsfd = (D0/(1+(p/pscale)+((p/pscale)**2))) # Updating the Single File Diffusion effects based on density at positions
        
        p_1[1:width] = p[1:width] +  dt * ( Dsfd[0:width-1]*p[0:width-1]/(dx**2) + Dsfd[2:width+1]*p[2:width+1]/(dx**2) - Dsfd[1:width]*2*p[1:width]/(dx**2)) ### WITH SFD 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] ### WITH SFD
                
        p_3[1:width] = p_2[1:width] +  dt * D0 * ( (p_2[0:width-1] + p_2[2:width+1] - 2*p_2[1:width]) / (dx**2))             ### NO SFD
        acetyl2[0:width] = acetyl2[0:width] + (acetylmultiplicity - acetyl2[0:width] ) * arate * dt * p_2[0:width] ### NO SFD
         
        telapsed += dt                     

        #counter2+=1

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        p[width]=p[width-1] # Closed Tube Boundary

        p_2,p_3 = p_3,p_2 # Updating concentration array ### FOR NO SFD CASE
        p_2[0] = p0 # resetting the boundary condition   ### FOR NO SFD CASE
        p_2[width]=p_2[width-1] # Closed Tube Boundary   ### FOR NO SFD CASE
       
        if (counter2%10==0):
                pTot[int(counter2/10)] = sum(p[1:width]) * dx
                aTot[int(counter2/10)] = sum(acetyl) * dx
                
                pTot2[int(counter2/10)] = sum(p_2[1:width]*dx) 
                aTot2[int(counter2/10)] = sum(acetyl2) * dx
                
                

                z = x/np.sqrt(4*D0*(telapsed))
                
                pTotAnalytical[int(counter2/10)] = p0*np.sqrt(4*D0*(telapsed)/np.pi)
                #pTotAnalytical[int(counter2/10)-1] = scii.quad(densityAnalytical,0,x[width-1],args=(telapsed),epsabs=1e-15)[0]
                #pTotAnalytical[int(counter2/10)-1] = sum(scis.erfc(z))
                acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
                aTotAnalytical[int(counter2/10)] = scii.quad(acetylAnalytical,0,x[width-1],args=(telapsed),epsabs=1e-15)[0] 
        counter2+=1
        #p_2[0] = p0
print(telapsed)
tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)

print('residual')
print(pTot2[0]-pTotAnalytical[0])
print('max residual')
print(max(abs(pTot2-pTotAnalytical)))

'''
### PLOTTING
parameterString = 'arate=%.2f \n L=%.2f \n dx=%.2f'%(arate,L,dx)

plt.scatter(x/L , p[0:width] , marker=',', s=2 , label='SFD') 
plt.scatter(x/L , p_2[0:width] , marker='*' , s=2 , label='No SFD')
plt.xlabel('x/L')
plt.ylabel('Density')
plt.legend()

plt.figure()
plt.scatter(x/L , acetyl , marker=',' , s=2 , label='SFD')
plt.scatter(x/L , acetyl2 , marker='*' , s=2 , label='No SFD')
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
