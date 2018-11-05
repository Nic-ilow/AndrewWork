import math as m
import numpy as np
import matplotlib.pyplot as plt
import random as ran
import argparse
import time
import pickle
import os
import datetime
import shutil

now = datetime.datetime.now()
dateFormat = 'Sim_{0}_{1}_{2}_{3}_{4}'.format(now.day,now.month,now.year,now.hour,now.minute)

if os.path.exists(dateFormat):
        shutil.rmtree(dateFormat)
os.makedirs(dateFormat)

parser=argparse.ArgumentParser(description='Command line inputs:') # needs to import argparse at top
parser.add_argument('-D0','--D0',default=2.7E5,help='diffusion as measured')	# Szyk measured diff
parser.add_argument('-dx','--dx',default=7.0,help='particle diam')		# nm (size of alphatTat1)
parser.add_argument('-koff','--koff',default=10.0,help='off rate (per second)')
parser.add_argument('-kon','--kon',default=1000.0,help='on rate (per second)')
parser.add_argument('-tmax','--tmax',default=32.0,help='maximum time in seconds')
parser.add_argument('-width','--width',default=16,help='width system in units of dx')
parser.add_argument('-arate','--arate',default=1.0,help='probability of particle acetylating site')
parser.add_argument('-am','--am',default=13,help='acetyl multiplicity, e.g. how many acetyl groups per site')
parser.add_argument('-tubules','--tubules',default=3,help='Number of simulations to save data for')
args=parser.parse_args()

D0=float(args.D0)
dx=float(args.dx)
koff=float(args.koff)
kon=float(args.kon)
tmax=float(args.tmax)
width=int(args.width)
arate=float(args.arate)
am = int(args.am)
tubuleSims = int(args.tubules)


if kon==0:
        fractionfree = 1.0
else:
        R = koff/kon
        fractionfree=R/(1.0+R)  # large R gives fractionfree=1 e.g. how often the particle is NOT bound
khop=2.0*D0/(dx*dx)/fractionfree 

# density 1 at 0
# density 0 at width (which isn't a site)

#
mt = np.arange(0,width)
mt2 = mt*dx
mtScaled = mt/(width-1.0)
global x,leftIn,rightOut,Nbound,Nfree,free,bound,x,asite,tElapsed

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
        global Nfree,tElapsed,asite,Nacetyl,Nbound,free,bound
        temp = free
        temp.extend(bound)
        pos = ran.choice(temp)
        temp2 = ran.random()
        if temp2<=(1-asite[pos]/am):
                asite[pos]+=1
                Nacetyl += 1
        tElapsed += dt

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

##############################################

tStart = time.time()
counter = 0

while counter<tubuleSims:
        print('Tubule_{0}'.format(counter))
        ran.seed(int(time.time()*1000000))
        x = np.zeros(width,int)  # 0 if empty, 1 if particle
        Density = np.zeros(np.size(x))

        asite = np.zeros(width,int) # 1 if acetylated, 0 if not
        Acetylation = np.zeros(np.size(asite))

        totPoints = []
        plotPoints = []
        NTot = []
        ATot = []
        pArray = []
        aArray = []

        Nbound = 0
        Nfree = 0
        Nacetyl = 0	# number of particles

        free = []
        bound = []

        rightOut = 0
        leftIn = 0

        tElapsed = 0

        plotCuts = [1, 2 , 4 , 8 , 16 , 30 ]
        netCuts = list(np.logspace(-200,0,num=200,base=1.1)*tmax)
        counter2 = 0
        while tElapsed <= tmax:
                fill() # Replacing the boundary

                totrate = Nfree*(kon+khop) + Nbound*koff +(Nbound+Nfree)*arate #Kinetic Monte Carlo
                dt = -1.0/totrate*m.log(1.0-ran.random())

                nextx = ran.random()*totrate

                if nextx < Nfree*kon:
                        binder('bind')

                elif nextx < Nfree*(kon+khop):
                        hop(ran.randrange(Nfree))

                elif nextx < Nfree*(kon+khop) + (Nbound+Nfree)*arate:
                        acetylate()

                else:
                        binder('unbind')

                Density += x*dt
                #Acetylation += asite*dt

                if any((tElapsed-t)>0 for t in plotCuts):
                        temp1 = np.copy(x)
                        temp2 = np.copy(asite)
                        pArray.append(temp1)
                        aArray.append(temp2)
                        plotPoints.append(plotCuts.pop(0))

                while any((tElapsed-t)>0 for t in netCuts):
                       NTot.append(sum(np.copy(x)))
                       ATot.append(sum(np.copy(asite)))
                       totPoints.append(netCuts.pop(0))

                counter2 += 1
        tEnd = time.time()
        tRun = (tEnd-tStart)/60
        print('Time Spent Running = %.2f'%tRun)
        print('Average Timestep = %2f'%(tmax/counter2))
        flux = [leftIn,rightOut]

        tubuleDens = 'DENS_{0}'.format(counter)
        tubuleAcetylation = 'ACETYL_{0}'.format(counter)
        tubuleNTot = 'NTOT_{0}'.format(counter)
        tubuleATot = 'ATOT_{0}'.format(counter)
        tubuleFlux = 'FLUX_{0}'.format(counter)
        np.save(os.path.join(dateFormat,tubuleDens),pArray)
        np.save(os.path.join(dateFormat,tubuleAcetylation),aArray)
        np.save(os.path.join(dateFormat,tubuleNTot),NTot)
        np.save(os.path.join(dateFormat,tubuleATot),ATot)
        np.save(os.path.join(dateFormat,tubuleFlux),flux)

        counter+=1

simInfo = [width,tmax,tubuleSims,arate,am,koff,kon]

np.savetxt(os.path.join(dateFormat,'SIMINFO'),simInfo,delimiter=',')
np.save(os.path.join(dateFormat,'TotTimes'),totPoints)
np.save(os.path.join(dateFormat,'SliceTimes'),plotPoints)

shutil.move(dateFormat,'Data')


