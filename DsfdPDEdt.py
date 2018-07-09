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
parser.add_argument('-width','--width',default=512)
parser.add_argument('-p0','--p0',default=1)
parser.add_argument('-dt','--dt',default=1E-6)
parser.add_argument('-arate','--arate',default=100)
parser.add_argument('-am','--am',default=1)
parser.add_argument('-pscale','--pscale',default=1000)
args = parser.parse_args()

global p0,arate,acetylmultiplicity,D0

### Relevant constants
D0 = float(args.D0)
dx = float(args.dx)
tmax = float(args.tmax)
width = int(args.width)
p0 = float(args.p0)
dt = float(args.dt)
arate = float(args.arate)
acetylmultiplicity = int(args.am)
pscale = float(args.pscale)

### Analytical Acetylation Function

def acetylAnalytical(x,telapsed):
        z = (x*dx) / np.sqrt(4*D0Scaled*telapsed)
        acetyl = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
        return  acetyl
        
#acetylAnalytical = np.vectorize(acetylAnalytical)


### Length and Steps
#L = 14336.0 #3 microns eg 3000 nanometers since we're in units of nm
L = 3000
#L = width*dx

### Our effective Diffusion Constant for taking derivative
D = D0/(dx**2)


### Plotting Cutoffs
base = 2
Samples = 5
tstart = time.time()
tmaxArray = np.logspace(-1,Samples,Samples+2,base=base)
#MaxDensityResidual = np.zeros( (np.size(dxArray),int(m.log(tmax,2)+2)) )

pArray = np.empty((np.size(tmaxArray),width))
pAnalyticArray = np.zeros_like(pArray)
pTotArray = np.zeros_like(pArray)
pTotAnalyticalArray = np.zeros_like(pArray)

DensityResidual = np.zeros( width )
AcetylResidual = np.zeros_like(DensityResidual)
NetDensityResidual = np.zeros_like(DensityResidual)
NetAcetylResidual = np.zeros_like(DensityResidual)

D0Scaled = (D0/(1+(p0/pscale)+((p0/pscale)**2)))
step=0
tTotstart = time.time()
D = D0/(dx**2)
#D = (D/(1+(p0/pscale)+((p0/pscale)**2))) #Single File Diffusion D
dt = (dx**2) / D0 * .5
print('dx=%.2f'%dx) 
print('D=%.3e'%D)
print('dt = %.3e'%dt)
x = np.arange(0,L,dx)
width = np.size(x)
t2 = 0
        
### Concentration Variables
p = np.zeros(width+1)
p_1 = np.zeros_like(p)

acetyl = np.zeros(width)
pTot = np.zeros(int(int(tmax/dt+1)/10)+1)
aTot = np.zeros_like(pTot)
pTotAnalytical = np.zeros_like(pTot)
aTotAnalytical = np.zeros_like(aTot)
### Boundary Condition
p[0] = p0
p[width-1] = 0

telapsed = 0
print(telapsed)
epsilon = dt/2

tstart=time.time()

counter = 0
counter2 = 0
tArray = np.arange(0,tmax,dt*10)
dt = (dx**2)/D0*.5
while telapsed<=tmax:
        Dsfd = (D/(1+(p/pscale)+((p/pscale)**2)))
        p_1[1:width] = p[1:width] + Dsfd[0:width-1] * dt * ( p[0:width-1] + p[2:width+1] - 2*p[1:width] ) 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

        
        telapsed += dt                     
        
        z = x/np.sqrt(4*D0Scaled*telapsed)
                            
        acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )                
        
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
                        
                        #densityExperiment[counter2,:] = p[0:width]
                        #acetylExperiment[counter2,:] = acetyl
                        print(telapsed)
                        counter2+=1

plt.figure(1)
plt.xlabel('x scaled (x/L)')
plt.ylabel('Density')
plt.legend()

plt.figure(2)
plt.xlabel('x scaled (x/L)')
plt.ylabel('Acetylation')
plt.legend()

DensityResidual = p[0:width] - scis.erfc(z)
AcetylResidual = acetyl-acetylationfit
NetDensityResidual =  pTot - pTotAnalytical
NetAcetylResidual = aTot - aTotAnalytical

tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)


tTotend=time.time()
trun = (tTotend-tTotstart)/60
print('Total Time Ran: %.2f'%trun)
parameterString = ' dt (s)=%.2e \np0=%.2f \nLength of Tubule (nm)=%.2f  '%(dt,p0,L)

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

