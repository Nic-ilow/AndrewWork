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
#L = 3E3 #3 microns eg 3000 nanometers since we're in units of nm
L = width*dx
x = np.arange(0,L,dx)
#width = np.size(x)

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

tstep = tmax/100000
tArray = np.arange(0,tmax,tstep)
xscaled = x/width/dx
pTot = np.empty(int(tmax/tstep)+1)
aTot = np.empty(int(tmax/tstep)+1)
analyticalpTot = np.empty_like(pTot)

counter = 0
tstart = time.time()
t1 = tstep
t2 = 0
t2Array = np.array( [0.5, 1, 2, 4, 8, 16, 32] )

pArray = np.empty((7,width))
pAnalyticArray = np.zeros_like(pArray)



Dsfd = (D/(1+(p0/pscale)+((p0/pscale)**2))) 
while telapsed<=tmax:

        p_1[1:width-1] = p[1:width-1] + D * ( p[0:width-2] + p[2:width] - 2*p[1:width-1] ) 
        acetyl[0:width] = acetyl[0:width] + (acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 

        
        telapsed += dt
        if abs(telapsed-(t2Array[t2]))<=epsilon:
                
                if t2==0:
                        z = x/np.sqrt(4*D0*telapsed)
                        a = 't= '+str(t2Array[0])+'s analytical'
                        b = 't= '+str(t2Array[0])+'s simulated'
                        c = 't= '+str(t2Array[0])+'s residual'
                        d = 'Acetylation t= '+str(t2Array[0])+'s residual'
                else: 
                        z = x/np.sqrt(4*D0*int(telapsed+.1))
                        a = 't= '+str(int(telapsed+.1))+'s analytical'
                        b = 't= '+str(int(telapsed+.1))+'s simulated'
                        c = 't= '+str(int(telapsed+.1))+'s residual'
                        d = 'Acetylation t= '+str(int(telapsed+.1))+'s residual'                       
                pArray[t2,:] = p_1[:]
                pAnalyticArray[t2,:] = scis.erfc(z)*p0
                 
                plt.figure(1)
                plt.plot(xscaled,pAnalyticArray[t2,:] , label=a) 
                plt.scatter(xscaled,pArray[t2,:],marker=',',s=2,label=b)
                
                residual = pAnalyticArray[t2,:]-pArray[t2,:]
                plt.figure(2)
                plt.plot(xscaled,residual,label=c)
                t2 +=1
                print(telapsed)    
                 
                acetylationfit = acetylmultiplicity * (1 - np.exp( -p0*arate*telapsed* ( (1+ (2 * (z**2) ) ) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) ))
                
                plt.figure(3)
                plt.scatter(xscaled,acetyl,label=b,marker=',',s=2)
                plt.plot(xscaled,acetylationfit,label=a)
                        
                       
                plt.figure(4)
                plt.scatter(xscaled,acetyl-acetylationfit,label=d,marker=',',s=2)

                        
        if abs(telapsed-t1)<=epsilon:
                pTot[counter] = sum(p)
                aTot[counter] = sum(acetyl)
                analyticalpTot[counter] = p0*np.sqrt(4*D0/(dx**2)*telapsed/np.pi)
                t1 += tstep
                counter += 1
                #print(telapsed)
        
        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        #p[width-1] = 0
        occupationDensity += p*dt
        acetylDensity += acetyl*dt
tend=time.time()
trun = (tend-tstart)/60 #In minutes
print(trun)
parameterString = 'dx (nm)=%.2f \ndt (s)=%.2e \nwidth (# of sites)=%.2f \np0=%.2f \nLength of Tubule (nm)=%.2f  '%(dx,dt,width,p0,L)
plt.figure(1)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Density')
plt.text(0.50,0.60, parameterString,verticalalignment='top',horizontalalignment='left')
plt.title('Density as a function of Position')
plt.legend()

plt.figure(2)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Density (p)')
#plt.text(0.05,0.95, parameterString,verticalalignment='top',horizontalalignment='left')
plt.title('Density Residual as a function of position')
plt.legend()

plt.figure(3)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Acetylation (A.U.)')
#plt.text(0.05,0.95, parameterString,verticalalignment='top',horizontalalignment='left')
plt.title('Acetylation as a function of Position')
plt.legend()

plt.figure(4)
plt.xlabel('Normalized Length (X/L)')
plt.ylabel('Acetylation (A.U.(')
#plt.text(0.05,0.95, parameterString,verticalalignment='top',horizontalalignment='left')
plt.title('Acetylation Residual as a function of Position')
plt.legend()

plt.figure(5)
plt.scatter(tArray,pTot,label='Total concentration',marker=',',s=2)
#plt.scatter(tArray,aTot,label='Total acetylation',marker=',',s=2)
plt.plot(tArray,analyticalpTot,label='Total Concentration Analytic')
plt.xlabel('time (s)')
plt.ylabel('Total Concentration')
#plt.text(0.05,0.95, parameterString,verticalalignment='top',horizontalalignment='left')
plt.title('Total concentration as a function of time')
plt.legend()

plt.figure(6)
plt.plot(tArray,pTot-analyticalpTot,label='total concentration residual')
plt.xlabel('time (s)')
plt.ylabel('Total Concentration residual')
plt.title('Total Concentration residual as a function of time')
#plt.text(0.05,0.95,parameterString,verticalalignment='top',horizontalalignment='left')
plt.show()


'''
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
'''
occupationDensity /= tmax*p0
acetylDensity /= tmax*acetylmultiplicity
'''
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
'''
