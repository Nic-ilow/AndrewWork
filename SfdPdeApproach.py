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
args = parser.parse_args()

### Relevant constants
width = int(args.width)
D0 = float(args.D0)
dx = float(args.dx)
tsteps = float(args.tsteps)
arate = float(args.arate)
p0 = float(args.p0)

p = np.zeros(width)
p_1 = np.zeros_like(p)
asite = np.zeros_like(p)
x = np.arange(width)
xscaled = x/(width-1)

### Density Variables
acetylDensity = np.zeros_like(p)
occupationDensity = np.zeros_like(p)


### Boundary Conditions
p[0]=p0

#J = -Dp *grad(p) = -D * p'
#dp/dt = -grad(J) = (D*p')'

telapsed = 0
dt = 1E-3
tmax = tsteps*dt
D = D0/(dx**2)
acetylMultiplicity = 5

#def secondOrder(p):
print ('dt:',dt, ' dx (in nm):',dx, ' width:',width, ' D0:',D0, ' D:',D)

def progress(p):
    newP1 = p[1] + D * ( p[0] - 2*p[1] + p[2] )
    return newP1


#np.vectorize(secondOrder)
np.vectorize(progress)

while telapsed<tmax:
    for i in range(width):
            
            randomAcetyl = rand.random()            

            if randomAcetyl <= arate*p[i]:
                    if asite[i] < acetylMultiplicity:
                            asite[i]+=1
                            #print 'what'
            if i==width-1:
                    pPass = np.array([p[width-2] , p[width-1] , 0])
            elif i==0:
                    pPass = np.array([p[0] , p[0] , p[1]])
            else:
                    pPass = np.array([p[i-1],p[i],p[i+1]])
            
            p_1[i] = progress(pPass)

    telapsed += dt
    p = p_1

    acetylDensity[:] += asite[:]*dt
    p[0] = p0
    p[width-1] = 0 
    occupationDensity[:] += p[:]*dt
occupationDensity /= tmax
acetylDensity /= tmax

normalizedAcetylDensity = acetylDensity / acetylMultiplicity
normalizedOccupationDensity = occupationDensity / p0
#plt.figure()
#plt.plot(p)

plt.figure()
plt.plot(normalizedAcetylDensity,color='r',label='Acetylation Density')
plt.plot(normalizedOccupationDensity,color='b',label='Occupation Density')
plt.ylabel('Normalized Time')
plt.xlabel('Normalized Length')
plt.legend()
plt.show()

