import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import pickle
import time
import math as m

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
D0 = float(args.D0)
dx = float(args.dx)
tmax = float(args.tmax)
width = int(args.width)
p0 = float(args.p0)
dt = float(args.dt)
arate = float(args.arate)
acetylmultiplicity = int(args.am)
pscale = float(args.pscale)

### Length and Steps
L = 14336.0 #3 microns eg 3000 nanometers since we're in units of nm
#L = width*dx
x = np.arange(0,L,dx)
width = np.size(x)

### Our effective Diffusion Constant for taking derivative
D = D0/(dx**2)
pscaleS = str(pscale)

### Plotting Cutoffs
epsilon = dt/2
t1 = 0

### Concentration Variables
p = np.zeros(width+1)
p_1 = np.zeros_like(p)


acetyl = np.zeros(width)
x=np.arange(width)*dx

### Density Tracking
occupationDensity = np.zeros(width)
acetylDensity = np.zeros(width)

### Boundary Condition
p[0] = p0
p[width-1] = 0

telapsed = 0

tArray = np.arange(0,tmax,dt)
xscaled = x/width/dx
pTot = np.empty(int(tmax/dt))
aTot = np.empty(int(tmax/dt))
analyticalpTot = np.empty_like(pTot)
testpTot = np.empty_like(pTot)

counter = 0
tstart = time.time()
t2 = 0
t2Array = np.geomspace(.5,tmax,m.log(tmax,2)+2)

pArray = np.empty((int((m.log(tmax,2)+2)),width))
pAnalyticArray = np.zeros_like(pArray)
MaxDensityResidual = np.zeros_like(t2Array)


Dsfd = (D/(1+(p0/pscale)+((p0/pscale)**2))) 
while telapsed<=tmax:

        p_1[1:width] = p[1:width] + D * dt * ( p[0:width-1] + p[2:width+1] - 2*p[1:width] ) 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

        
        telapsed += dt
        if (abs(telapsed-(2**(t2-1)))<=epsilon):
        
                z = x/np.sqrt(4*D0*int(telapsed+.1))
                a = 't= '+str(int(telapsed+.1))+'s analytical'
                b = 't= '+str(int(telapsed+.1))+'s simulated'
                c = 't= '+str(int(telapsed+.1))+'s residual'
                d = 'Acetylation t= '+str(int(telapsed+.1))+'s residual'                       
                pArray[t2,:] = p_1[0:width]
                pAnalyticArray[t2,:] = scis.erfc(z)*p0
                 
                ### PLOTTING DENSITY  
                plt.figure(1)
                plt.plot(xscaled,pAnalyticArray[t2,:] , label=a) 
                plt.scatter(xscaled,pArray[t2,:],marker=',',s=2,label=b)
                
                residual = pAnalyticArray[t2,:]-pArray[t2,:]
                MaxDensityResidual[t2] = max(abs(residual))
                
                ### PLOTTING RESIDUAL
                plt.figure(2)
                plt.plot(xscaled,residual,label=c)
                
                t2 +=1
                print(telapsed)    
                 
                acetylationfit = acetylmultiplicity * (1 - np.exp( -p0*arate*telapsed* ( (1+ (2 * (z**2) ) ) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ))
                
                ### PLOTTING ACETLYATION
                plt.figure(3)
                plt.scatter(xscaled,acetyl,label=b,marker=',',s=2)
                plt.plot(xscaled,acetylationfit,label=a)
                        
                ### PLOTTING RESIDUAL
                plt.figure(4)
                plt.scatter(xscaled,acetyl-acetylationfit,label=d,marker=',',s=2)
                print(telapsed)

                        
        y=np.sqrt(4*D0*telapsed)
        pTot[counter] = sum(p_1[0:width])
        aTot[counter] = sum(acetyl)
        analyticalpTot[counter] = p0*np.sqrt(4*D*telapsed/np.pi)
        
        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        p[width]=p[width-1] # Making it a closed end
        #p[width-1] = 0
        occupationDensity += p[0:width]*dt
        acetylDensity += acetyl*dt
        counter+=1
tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)

parameterString = 'dx (nm)=%.2f \ndt (s)=%.2e \nwidth (# of sites)=%.2f \np0=%.2f \nLength of Tubule (nm)=%.2e  '%(dx,dt,width,p0,L)
'''
plt.figure(1)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Density')
plt.text(0.50,0.60, parameterString,verticalalignment='top',horizontalalignment='left')
plt.title('p(x,t)')
plt.legend()

plt.figure(2)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Density (p)')
plt.title('p(x,t) Residual')
plt.legend()

plt.figure(3)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Acetylation (A.U.)')
plt.title('a(x,t)')
plt.legend()

plt.figure(4)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Acetylation (A.U.(')
plt.title('a(x,t Residual)')
plt.legend()
'''
plt.figure(5)
plt.loglog(tArray,pTot,label='Total concentration',marker=',')
#plt.scatter(tArray,aTot,label='Total acetylation',marker=',',s=2)
plt.loglog(tArray,analyticalpTot,label='Total Concentration Analytic')
#plt.plot(tArray,testpTot,label='Total Concentration Analytic')
plt.xlabel('time (s)')
plt.ylabel('Total Concentration')
plt.title('N(t)')
plt.legend()

plt.figure(6)
plt.loglog(tArray,pTot-analyticalpTot,label='total concentration residual')
plt.xlabel('time (s)')
plt.ylabel('Total Concentration residual')
plt.title('N(t) residual')
plt.show()
'''
#residualString = ('MaxDensityResidualDxis%.1i'%dx)
residualString = ('MaxDensityResidualDtis%.1e'%dt)
np.save(residualString,MaxDensityResidual)
'''
