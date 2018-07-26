import math as m
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import scipy.integrate as scii
import pickle
import time

parser = argparse.ArgumentParser(description='Command line inputs:')
#parser.add_argument('-D','--D',default=2.7E-13)
#parser.add_argument('-dx','--dx',default=7.0)
#parser.add_argument('-tmax','--tmax',default=1.0)
#parser.add_argument('-p0','--p0',default=1)
#parser.add_argument('-arate','--arate',default=100)
#parser.add_argument('-am','--am',default=1)
#parser.add_argument('-pscale','--pscale',default=1000)
#parser.add_argument('-dt','--dt',default=1e-6)

#parser.add_argument('-du','--du',default=.04) ### SCALED OUT LENGTH STEP
#parser.add_argument('-uLength','--uLength',default=17) ### SCALED OUT LENGTH OF SYSTEM
#parser.add_argument('-dv','--dv',default=8e-7) ### SCALED OUT TIME STEP
#parser.add_argument('-vmax','--vmax',default=0.8) ### SCALED OUT TMAX
#parser.add_argument('-pscale','--pscale',default=1000.0)
#parser.add_argument('-p0','--p0',default=1.0)
#parser.add_argument('-dx','--dx',default=1e-9)
parser.add_argument('-p0hat','--p0hat',default=1e-4)
parser.add_argument('-dxbar','--dxbar',default=6.0)
parser.add_argument('-xbarmax','--xbarmax',default=1800.0)
parser.add_argument('-dtbar','--dtbar',default=.1)
parser.add_argument('-tbarmax','--tbarmax',default=1.0e5)
args = parser.parse_args()

p0hat = float(args.p0hat)
dxbar = float(args.dxbar)
xbarmax = float(args.xbarmax)
dtbar = float(args.dtbar)
tbarmax = float(args.tbarmax)

xbar = np.arange(0,xbarmax,dxbar)
tbar = np.arange(0,tbarmax,dtbar)

epsilon = dtbar/2
width = np.size(xbar)       

### pTot and aTot plotting array
tArray = np.arange(1*dtbar,tbarmax+dtbar,dtbar*10)

### Initializing counters and a timer 
tstart = time.time()
counter=0
### Array Pre-allocation
p = np.zeros(width+1)
p_1 = np.zeros_like(p)
telapsed = 0

acetyl = np.zeros(width)
pTot = np.zeros((int(int(tbarmax/dtbar+1)/10+1)))
aTot = np.zeros_like(pTot)
aTotAnalytical =  np.zeros_like(aTot)
pTotAnalytical = np.zeros_like(aTotAnalytical)

### Boundary Condition
p[0] = p0hat
p[width-1] = 0
counter2 = 0

while telapsed<=tbarmax: # Iterating the system of tmax amount of seconds
        pscaler = 1+p+p**2     
        p_1[1:width] = p[1:width] +   dtbar *  ( (p[0:width-1]/pscaler[0:width-1] + p[2:width+1]/pscaler[2:width+1] - 2*p[1:width]/pscaler[1:width]) / (dxbar**2))  
        acetyl[0:width] = acetyl[0:width] + (1 - acetyl[0:width] ) * dtbar * p[0:width] ### NO SFD
         
        telapsed += dtbar                     

        #counter2+=1

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0hat # resetting the boundary condition
        p[width]=p[width-1] # Closed Tube Boundary

        if (counter2%10==0):
                pTot[int(counter2/10)] = sum(p[1:width]) * dxbar
                aTot[int(counter2/10)] = sum(acetyl) * dxbar
                #z = xbar/np.sqrt(4*D*(telapsed))
                
                #pTotAnalytical[int(counter2/10)] = p0*np.sqrt(4*D*(telapsed)/np.pi)
                #pTotAnalytical[int(counter2/10)-1] = scii.quad(densityAnalytical,0,x[width-1],args=(telapsed),epsabs=1e-15)[0]
                #pTotAnalytical[int(counter2/10)-1] = sum(scis.erfc(z))
                #acetylationfit = 1 - np.exp( -p0 * arate * telapsed *( (1+2*(z**2)) * scis.erfc(z) - 2*z*np.exp(-(z**2))/np.sqrt(np.pi) ) )
                #aTotAnalytical[int(counter2/10)] = scii.quad(acetylAnalytical,0,xbar[width-1],args=(telapsed),epsabs=1e-15)[0] 
        counter2+=1
        #p_2[0] = p0
print((time.time()-tstart)/60)
#counter+=1
#print(telapsed)
#tend=time.time()
#trun = (tend-tstart)/60 #In minutes
#print(trun)

#print('residual')
#print(pTot2[0]-pTotAnalytical[0])
#print('max residual')
#print(max(abs(pTot2-pTotAnalytical)))
plt.scatter(xbar[1:width],p[1:width],s=1)
plt.show()
'''
### PLOTTING
parameterString = 'arate=%.2f \n L=%.2f \n dx=%.2f'%(arate,L,dx)
pTotResidual = np.zeros_like(dxArray)
#
plt.scatter(x/L , p[0:width] , marker=',', s=2 , label='SFD') 
plt.scatter(x/L , p_2[0:width] , marker='*' , s=2 , label='No SFD')
plt.xlabel('x/L')
plt.ylabel('Density')
plt.legend()

plt.figure()
plt.scatter(x/L , acetyl , marker=',' , s=2 , label='SFD')
plt.scatter(x/L , acetyl , marker='*' , s=2 , label='No SFD')
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
