import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.special as scis
import pickle

parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-D0','--D0',default=2.7E5)
parser.add_argument('-dx','--dx',default=7.0)
parser.add_argument('-tmax','--tmax',default=2.0)
parser.add_argument('-width','--width',default=32)
parser.add_argument('-p0','--p0',default=1)
parser.add_argument('-dt','--dt',default=1E-6)
parser.add_argument('-arate','--arate',default=100)
parser.add_argument('-am','--am',default=1)
parser.add_argument('-pscale','--pscale',default=100)
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
D = 2*D0*dt/(dx**2)
pscaleS = str(pscale)
### Plotting Cutoffs
epsilon = dt/1000
t0 = 1
t1 = 2
t2 = 4
t3 = 8
t4 = 16

### Concentration Variables
#p = np.random.rand(width)
p = np.zeros(width)
ptest = p
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

### Not Needed Anymore
'''
def progress(p):

        newP1 = p[1] + ( D ( p[0] - 2*p[1] + p[2] ) )
        return newP1

def rhopass(i):
        if i==width-1:
                pPass = np.array( [ p[i-1] , p[i] , p[i] ] )
        elif i==0:
                pPass = np.array( [ p[i] , p[i] , p[i+1] ] )
        else:
                pPass = np.array( [ p[i-1] , p[i] , p[i+1] ] )
        
        return pPass
np.vectorize(progress)
'''

x2 = x/width/dx
while telapsed<=tmax:

        p_1[1:width-1] = p[1:width-1] + (D/(1+(p0/pscale)+((p0/pscale)**2))) * ( p[0:width-2] + p[2:width] - 2*p[1:width-1] ) 
        acetyl[0:width] = acetyl[0:width] + ( acetylmultiplicity - acetyl[0:width] ) * arate * dt * p[0:width] 
        telapsed += dt
        
        '''
        if abs(telapsed-t0)<epsilon:
                #plt.figure()
                plt.plot(x2,occupationDensity/telapsed,label=('t= t0')) 
        #        plt.legend()
                print(telapsed)
        if abs(telapsed-t1)<epsilon:
                #plt.figure()
                plt.plot(x2,occupationDensity/telapsed,label=('t= t1'))
         #       plt.legend()
                print(telapsed)
        if abs(telapsed-t2)<epsilon:
                #plt.figure()
                plt.plot(x2,occupationDensity/telapsed,label=('t= t2'))
          #      plt.legend()
                print(telapsed)
        if abs(telapsed-t3)<epsilon:
                #plt.figure()
                plt.plot(x2,occupationDensity/telapsed,label=('t= t3'))
           #     plt.legend()
                print(telapsed)
        if abs(telapsed-t4)<epsilon:
                #plt.figure()
                plt.plot(x2,occupationDensity/telapsed,label=('t= t4'))
            #    plt.legend()
                print(telapsed)
        '''

        p,p_1 = p_1,p # Updating concentration array
        p[0] = p0 # resetting the boundary condition
        #p[width-1] = 0
        occupationDensity += p*dt
        acetylDensity += acetyl*dt
test2 = np.linspace(p0,0,width)
test = scis.erfc(x/np.sqrt(4*D0*tmax))*test2/p0

occupationDensity /= tmax*p0
acetylDensity /= tmax
#np.save('acetyl'+pscaleS,acetylDensity)
#np.save('occupation'+pscaleS,occupationDensity)

plt.figure()
plt.scatter(x2,occupationDensity-test,label='residual')
plt.figure()
plt.plot(x2,occupationDensity,label='t= tmax')
plt.plot(x2,test,label='tmax erfc fit')
plt.plot(x2,acetylDensity,label='acetyl occupancy')
plt.legend()
plt.show()

