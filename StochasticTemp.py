import math as m
import numpy as np
import matplotlib.pyplot as plt
import random as ran
import argparse
import time
import pickle

parser=argparse.ArgumentParser(description='Command line inputs:') # needs to import argparse at top
parser.add_argument('-R','--R',default=0.01,help='koff/kon ratio')
parser.add_argument('-D0','--D0',default=2.7E5,help='diffusion as measured')	# Szyk measured diff
parser.add_argument('-dx','--dx',default=7.0,help='particle diam')		# nm (size of alphatTat1)
parser.add_argument('-koff','--koff',default=10.0,help='off rate (per second)')
parser.add_argument('-tmax','--tmax',default=32.0,help='maximum time in seconds')
parser.add_argument('-nrun','--nrun',default=0,help='number of this run')
parser.add_argument('-name','--name',default='MT',help='datafile name base in single quotes')
parser.add_argument('-seed','--seed',default=0,help='seed for random number')
parser.add_argument('-width','--width',default=16,help='width system in units of dx')
parser.add_argument('-arate','--arate',default=1,help='probability of particle acetylating site')
parser.add_argument('-am','--am',default=13,help='acetyl multiplicity, e.g. how many acetyl groups per site')
args=parser.parse_args()

R=float(args.R)
D0=float(args.D0)
dx=float(args.dx)
koff=float(args.koff)
tmax=float(args.tmax)
nrun=int(args.nrun)
name=str(args.name)
seed=int(args.seed)+int(time.time()*1000000)   # so we essentially use the time, but we allow external offset
width=int(args.width)
arate=float(args.arate)
am = int(args.am)

ran.seed(seed)
fractionfree=R/(1.0+R)  # large R gives fractionfree=1 e.g. how often the particle is NOT bound
khop=2.0*D0/(dx*dx)/fractionfree 

kon = koff/R  # (per second) 0.0 gives SD (singular limit)  (0.001 gives 10/s)  R=koff/kon

# density 1 at 0
# density 0 at width (which isn't a site)

print("R:",R," D0:",D0," dx:",dx," kon:",kon," koff:",koff," khop:", khop)
print("tmax:",tmax," nrun:",nrun," name:",name," seed:",seed," width:",width," arate:",arate)

#
mt = np.arange(0,width)
mtScaled = mt/(width-1.0)
global x,leftIn,rightOut,Nbound,Nfree,free,bound,x,asite,tElapsed
x = np.zeros(width,int)  # 0 if empty, 1 if particle
Density = np.zeros(np.size(x))

asite = np.zeros(width,int) # 1 if acetylated, 0 if not
Acetylation = np.zeros(np.size(asite))

pArray = []
aArray = []

Nbound = 0
Nfree = 0
Nacetyl = 0	# number of particles
free = []
bound = []

rightOut = 0
leftIn = 0

###################################################
def hop(this):
        global x,free,rightOut,Nfree,tElapsed
        pos = free[this]

        if ran.random()<0.5:
                if pos>0 and x[pos-1]==0:
                        x[pos] = 0
                        x[pos-1] = 1
                        free[this] = pos-1
        else:
                if pos==width-1:
                        rightOut += 1
                        x[pos] = 0
                        free.pop(this)
                        Nfree-=1
                elif x[pos+1] == 0:
                        x[pos] = 0
                        x[pos+1] = 1
                        free[this] = pos+1
        tElapsed += dt
##############################################
def acetylate():
        global Nfree,tElapsed,asite,Nacetyl,Nbound
        temp = ran.randrange(Nfree+Nbound)
        if temp<Nfree:
                pos = free[temp]
                if asite[pos]<am:
                        asite[pos]+=1
                        tElapsed += dt 
                        Nacetyl += 1
        else:
                pos = bound[temp-Nfree]
                if asite[pos]<am:
                        asite[pos]+=1
                        tElapsed += dt
                        Nacetyl += 1
###############################################
def binder(discriminator):
        global free,bound,Nbound,Nfree,tElapsed
        if discriminator=='unbind':
                pos = bound.pop(ran.randrange(Nbound))
                free.append(pos)
                Nbound -= 1 
                Nfree += 1
        
        else:
                pos = free.pop(ran.randrange(Nfree))
                bound.append(pos)
                Nfree -= 1
                Nbound += 1

        tElapsed += dt

###############################################
def fill():
        global Nfree,leftIn,free,x
        if x[0] == 0:
                x[0] = 1
                free.append(0)
                Nfree += 1
                leftIn +=1



tElapsed = 0
plotCuts = np.array( [tmax/4 , tmax/2 , 3*tmax/4 , tmax] )
plotPoints = []
tStart = time.time()
while tElapsed < tmax:
        
        fill()
        tTemp = tElapsed
        totrate = Nfree*(kon+khop) + Nbound*koff +(Nbound+Nfree)*arate
        dt = -1.0/totrate*m.log(1.0-ran.random())
        
        eps = dt/2

        nextx = ran.random()*totrate

        if nextx < Nfree*kon:
                binder('bind')
        elif nextx < Nfree*(kon+khop):
                hop(ran.randrange(Nfree))
        
        elif nextx < Nfree*(kon+khop) + (Nbound+Nfree)*arate:
                acetylate()
        
        else:
                binder('unbind')
        
        if tTemp!=tElapsed:
                Density += x*dt
                Acetylation += asite*dt

        
        if any(abs(tElapsed-plotCuts)<=eps):
                aArray.append(Acetylation)
                pArray.append(Density)
                plotPoints.append(plotCuts.pop(0))
                print('saved at time = %.2f'%tElapsed)
tEnd = time.time()
tRun = (tEnd-tStart)/60

print('Time Spent Running = %.2f'%tRun)
plotPoints = np.array(plotPoints)
for i , value in enumerate(plotCuts):
        plt.figure(1)
        plt.plot(mtScaled,pArray[i]/plotPoints[i],label=('t = %.2f s'%value))
        
        plt.figure(2)
        plt.plot(mtScaled,aArray[i]/am/plotPoints[i],label=('t = %.2f s'%value))

plt.figure(1)
plt.xlabel('Scaled Length')
plt.ylabel('Density')
plt.legend()

plt.figure(2)
plt.xlabel('Scaled Length')
plt.ylabel('Acetylation')
plt.legend()

plt.show()
