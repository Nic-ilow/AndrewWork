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
        z = (x*dx) / np.sqrt(4*D0*telapsed)
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
tstart = time.time()

pArray = np.empty((7,width))
pAnalyticArray = np.zeros_like(pArray)
dxArray = np.arange(1,1.1,.1)
#MaxDensityResidual = np.zeros( (np.size(dxArray),int(m.log(tmax,2)+2)) )

MaxDensityResidual = np.zeros( np.size(dxArray) )
MaxAcetylResidual = np.zeros_like(MaxDensityResidual)
MaxNetDensityResidual = np.zeros_like(MaxDensityResidual)
MaxNetAcetylResidual = np.zeros_like(MaxDensityResidual)

Dsfd = (D/(1+(p0/pscale)+((p0/pscale)**2))) 
step=0
tTotstart = time.time()
for dx in dxArray:
        D = D0/(dx**2)
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
        while telapsed<=tmax:

                if abs(telapsed-(tmax/4))<=epsilon or abs(telapsed-tmax/2)<=epsilon or abs(telapsed-tmax/10) <= epsilon:
                        print(telapsed)

                p_1[1:width] = p[1:width] + D * dt * ( p[0:width-1] + p[2:width+1] - 2*p[1:width] ) 
                acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

                
                telapsed += dt                     
                
                z = x/np.sqrt(4*D0*telapsed)
                                    
                acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )                
                
                p,p_1 = p_1,p # Updating concentration array
                p[0] = p0 # resetting the boundary condition
                p[width]=p[width-1] # Closed Tube Boundary
                
                if (counter%10)==0:
                        pTot[int(counter/10)] = sum(p[0:width])
                        pTotAnalytical[int(counter/10)] = p0 * np.sqrt(4*D*telapsed/np.pi)
                
                        aTot[int(counter/10)] = sum(acetyl)
                        aTotAnalytical[int(counter/10)] = scii.quad(acetylAnalytical,0,L,args=(telapsed))[0] 
                counter+=1
        
        acetylResidual = acetylationfit - acetyl
        residual = scis.erfc(z)*p0-p_1[0:width]
        pTotResidual = pTot-pTotAnalytical
        aTotResidual = aTot-aTotAnalytical
        
        MaxDensityResidual[step] = max(abs(residual))
        MaxAcetylResidual[step] = max(abs(acetylResidual))
        MaxNetDensityResidual[step] = max(abs(pTotResidual))
        MaxNetAcetylResidual[step] = max(abs(aTotResidual))
               
        step += 1          
        tend=time.time()
        trun = (tend-tstart)/60 #In minutes
        print(trun)


tTotend=time.time()
trun = (tTotend-tTotstart)/60
print('Total Time Ran: %.2f'%trun)
parameterString = ' dt (s)=%.2e \np0=%.2f \nLength of Tubule (nm)=%.2f  '%(dt,p0,L)

plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dxArray,MaxDensityResidual,label=('Residual at t=%.2f'%tmax))
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Density')
plt.legend()


plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dxArray,MaxAcetylResidual,label=('Residual at t=%.2f'%tmax))
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Acetylation')
plt.legend()

fig = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dxArray,MaxNetDensityResidual)
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Total Density')

plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dxArray,MaxNetAcetylResidual)
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Total Acetylation')
plt.show()
'''
plt.figure()
plt.loglog(dxArray,MaxDensityResidual[:,1],label='Residual at t=1.0s')
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Density')
plt.legend()

plt.figure()
plt.loglog(dxArray,MaxDensityResidual[:,2],label='Residual at t=2.0s')
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Density')
plt.legend()

plt.figure()
plt.title(parameterString)
plt.loglog(dxArray,MaxDensityResidual[:,3],label='Residual at t=4.0s')
plt.xlabel('dx (nm)')
plt.ylabel('Max Residual of Density')
plt.legend()
plt.show()
'''
