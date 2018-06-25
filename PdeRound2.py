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
epsilon = dt/2
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
tplot=1
steadyStateP = np.linspace(p0,0,width)

t2 = 0
t2Array = np.array( [1, 2, 4, 8, 16, 32] )

pArray = np.empty_like(t2Array)
pAnalyticArray = np.empty_like(t2Array)



Dsfd = (D/(1+(p0/pscale)+((p0/pscale)**2))) 
while telapsed<=tmax:

        p_1[1:width-1] = p[1:width-1] + D * ( p[0:width-2] + p[2:width] - 2*p[1:width-1] ) 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

        
        if abs(telapsed-t1)<epsilon:
                if abs(telapsed-tplot)<=epsilon:
                        z = x/np.sqrt(4*D0*telapsed)
                        if abs(telapsed-(t2Array[t2]))<=epsilon:
                                pArray[t2] = p
                                pAnalyticArray[t2] = scis.erfc(z) 
                                t2 +=1
                                print(telapsed)
                        
                        a = 't= '+str(tplot)+'s analytical'
                        b = 't= '+str(tplot)+'s simulated'
                        c = 'Occupation t= '+str(tplot)+'s residual'
                        d = 'Acetylation t= '+str(tplot)+'s residual'
                        
                        
                        plt.figure(1)
                        plt.scatter(xscaled,p,label=b,marker=',',s=2)
                        plt.plot(xscaled,scis.erfc(z),label=a)


                        plt.figure(2) 
                        acetylationfit = acetylmultiplicity * (1 - np.exp( -steadyStateP*arate*telapsed* ( (1+ (2 * (z**2) ) ) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ))
                        plt.scatter(xscaled,acetyl,label=b,marker=',',s=2)
                        plt.plot(xscaled,acetylationfit,label=a)
                        
                        
                        plt.figure(3)
                        plt.scatter(xscaled,p-(scis.erfc(z)*steadyStateP),label=c,marker=',',s=2)
                        
                        
                        plt.figure(4)
                        plt.scatter(xscaled,acetyl-acetylationfit,label=d,marker=',',s=2)
                        
                        
                        tplot += 1
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

plt.figure(1)
plt.xlabel('Normalized Length of Microtubule')
plt.ylabel('Concentration')
plt.legend()
plt.figure(2)
plt.xlabel('Normalized Length of Microtubule')
plt.ylabel('Acetylation')
plt.legend()
plt.figure(3)
plt.xlabel('Normalized Length of Microtubule')
plt.ylabel('Concentration')
plt.legend()
plt.figure(4)
plt.xlabel('Normalized Length of Microtubule')
plt.ylabel('Acetylation')
plt.legend()
occupationDensity /= tmax*p0
acetylDensity /= tmax*acetylmultiplicity

#np.save('acetyl'+pscaleS,acetylDensity)
#np.save('occupation'+pscaleS,occupationDensity)
z = x/np.sqrt(4*D0*tmax)

occupationFit = scis.erfc(z) #* steadyStateP
acetylationfit = acetylmultiplicity * (1 - np.exp( -arate*tmax* ( (1+ (2 * (z**2) ) ) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ))
plt.legend()

plt.figure()
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
plt.scatter(tArray,pTot,label='Total concentration',marker=',',s=2)
plt.scatter(tArray,aTot,label='Total acetylation',marker=',',s=2)
plt.xlabel('time (s)')
plt.ylabel('Total Concentration/Acetlyation throughout tubule')
plt.legend()
plt.show()

