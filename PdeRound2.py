import Taylor as Bunny
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
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

### Relevant constants
width = int(args.width)
D0 = float(args.D0)
dx = float(args.dx)
tmax = float(args.tmax)
p0 = float(args.p0)
dt = float(args.dt)
arate = float(args.arate)
acetylmultiplicity = int(args.am)
pscale = float(args.pscale)

### Our effective Diffusion Constant for taking derivative
D = D0*dt/(dx**2)
pscaleS = str(pscale)

### Plotting Cutoffs
epsilon = dt/10
t1 = 0

### Concentration Variables
p = np.zeros(width)
p_1 = np.zeros_like(p)


acetyl = np.zeros(width)
x=np.arange(width)*dx

### Density Tracking
occupationDensity = np.zeros_like(p)
acetylDensity = np.zeros_like(p)

### Boundary Condition
p[0] = p0
p[width-1] = 0

telapsed = 0

tstep = tmax/10000
tArray = np.arange(0,tmax+tstep,tstep)
xscaled = x/width/dx
pTot = np.empty(int(tmax/tstep)+1)
aTot = np.empty(int(tmax/tstep)+1)
counter = 0
tstart = time.time()
while telapsed<=tmax:

        p_1[1:width-1] = p[1:width-1] + (D/(1+(p0/pscale)+((p0/pscale)**2))) * ( p[0:width-2] + p[2:width] - 2*p[1:width-1] ) 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

        
        if abs(telapsed-t1)<epsilon:
                #if t1!=-0:
                        #plt.plot(xscaled,occupationDensity/telapsed,label=('t= t1'))
                pTot[int(counter/tstep)] = sum(p)
                aTot[int(counter/tstep)] = sum(acetyl)
                t1 += tstep
                counter += tstep
                #print(telapsed)
        
        telapsed += dt
        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        #p[width-1] = 0
        occupationDensity += p*dt
        acetylDensity += acetyl*dt
tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)

steadyStateP = np.linspace(p0,0,width)
occupationDensity /= tmax*p0
acetylDensity /= tmax*acetylmultiplicity

#np.save('acetyl'+pscaleS,acetylDensity)
#np.save('occupation'+pscaleS,occupationDensity)
z = x/np.sqrt(4*D0*tmax)

occupationFit = scis.erfc(z) * steadyStateP
acetylationfit = acetylmultiplicity * (1 - np.exp( -steadyStateP*arate*tmax* ( (1+ (2 * (z**2) ) ) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ))

plt.plot(xscaled,occupationDensity,label='occupation average density')
plt.plot(xscaled,acetylDensity,label='acetylation average density')
plt.legend()

plt.figure()
plt.plot(xscaled,p,label='occupation at tmax')
plt.plot(xscaled,occupationFit,label='occupation fit at tmax')
plt.legend()

plt.figure()
plt.plot(xscaled,(p-occupationFit),label='occupation residual at tmax')
plt.plot(xscaled,acetyl-acetylationfit,label='acetyl residual at tmax')
plt.legend()

plt.figure()
plt.plot(xscaled,acetyl,label='acetylation at tmax')
plt.plot(xscaled,acetylationfit,label='acetyl fit at tmax')
plt.legend()

plt.figure()
plt.plot(tArray,pTot,label='Total rho')
plt.plot(tArray,aTot,label='Total acetylated')
plt.legend()
plt.show()

