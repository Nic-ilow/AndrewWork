import os
import time
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time as t
import argparse
import random as rand
parser = argparse.ArgumentParser(description='Command line inputs:')
parser.add_argument('-D0','--D0',default=.27,help='diffusion as measured')
parser.add_argument('-dx','--dx',default=7.0,help='particle diam')
parser.add_argument('-tsteps','--tsteps',default=1000,help='maximum time in seconds')
parser.add_argument('-width','--width',default=16,help='width system in units of dx')
parser.add_argument('-arate','--arate',default=100,help='probability of acetylation when a particle is present per second')
parser.add_argument('-p0','--p0',default=1,help='Source concentration')
parser.add_argument('-kon','--kon',default=1000,help='Bind Rate per second for one particle')
parser.add_argument('-koff','--koff',default=10,help='Unbind rate per second for one particle')
parser.add_argument('-am','--am',default=1,help='acetylation Multiplicity')
args = parser.parse_args()

### Relevant constants
width = int(args.width)
D0 = float(args.D0)
dx = float(args.dx)
tsteps = float(args.tsteps)
arate = float(args.arate)
p0 = float(args.p0)
koff = float(args.koff)
kon = float(args.kon)
acetylMultiplicity = int(args.am)
seed = int(time.time()*1000000)
rand.seed(seed)
global bound
global asite

### Relevant Variables
p = np.zeros(width)    #Old Concentration
p_1 = np.zeros_like(p) #Updated Concentration
asite = np.zeros_like(p) #0 means no Acetylated , numbers above means number of acetylated sites
x = np.arange(width)
xscaled = x/(width-1)

### Density Variables
acetylDensity = np.zeros_like(p)
occupationDensity = np.zeros_like(p)
boundDensity = np.zeros_like(p)

### Bound / Unbound state tracking
bound = np.zeros(width+1)
bound[width] = 1


### Boundary Conditions
p[0]=p0

#J = -Dp *grad(p) = -D * p'
#dp/dt = -grad(J) = (D*p')'

telapsed = 0
dt = 1E-6
tmax = tsteps*dt
D = D0/(dx**2)

print ('dt:',dt, ' dx (in nm):',dx, ' width:',width, ' D0:',D0, ' D:',D)


def progress(p):
   
        newP1 = p[1] + D * ( p[0] - 2*p[1] + p[2] )
        return newP1


def rhoPass(i):
        
        if bound[i]==1 or (bound[i-1]==1 and bound[i+1]==1):
                pPass = np.array( [ p[i] , p[i] , p[i] ] )

        elif bound[i-1]==1:
                pPass = np.array( [ p[i] , p[i] , p[i+1] ] )
        
        elif bound[i+1]==1:
                pPass = np.array( [ p[i-1] , p[i] , p[i] ] )
        
        else:
                pPass = np.array( [ p[i-1] , p[i] , p[i+1] ] )

        return pPass

bindcount = 0
unbindcount=0
def binder(bound): #1 means Bound, 0 means unbound
        global unbindcount
        global bindcount
        for j in range(width):
                if ( bound[j]==0 and (rand.random() <= kon*dt) ): #If unbound, and being bound
                        bound[j]=1
                        bindcount +=1
                
                elif ( bound[j]==1 and (rand.random() <= koff*dt)) : #If bound, and being unbound
                        bound[j]=0
                        unbindcount +=1
# Note , bindcount - unbindcount should be sum(Bound) - 1
        return bound 

def acetylate(pos): #Checks if the position will be acetylated at a given timestep 

        if asite[pos] < acetylMultiplicity:
                if rand.random() <= arate*p[pos]*dt:
                        asite[pos]+=1

#np.vectorize(secondOrder)
np.vectorize(progress)
np.vectorize(binder)

while telapsed<tmax: #Iterate over some total time

        for i in range(width):
                acetylate(i)
                bound = binder(bound)
                if bound[i] != 1: #No point in doing computations if site is bound, becase p won't change there
                        p_1[i] = progress(rhoPass(i))


    
        telapsed += dt
        p = p_1
    
        boundDensity += bound[0:width]*dt
        acetylDensity += asite*dt
        p[0] = p0 
        occupationDensity += p*dt

occupationDensity /= tmax
acetylDensity /= tmax
boundDensity /= tmax

normalizedAcetylDensity = acetylDensity / acetylMultiplicity
normalizedOccupationDensity = occupationDensity / p0
#plt.figure()
#plt.plot(p)
os.system('say "yo shit done"')
plt.figure()
plt.plot(normalizedAcetylDensity,color='r',label='Acetylation Density')
plt.plot(normalizedOccupationDensity,color='b',label='Occupation Density')
plt.plot(boundDensity,color='g',label='Time Spent Bound')
plt.ylabel('Normalized Time')
plt.xlabel('Normalized Length')
plt.legend()
plt.title('test')
plt.show()
