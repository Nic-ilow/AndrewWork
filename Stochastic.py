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
parser.add_argument('-arate','--arate',default=100,help='probability of particle acetylating site')
paser.add_argument('-am','--am',default=1,help='acetyl multiplicity, e.g. how many acetyl groups per site')
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
khop=2.0*D0/(dx*dx)/fractionfree  # 2e6, correction factor is to use bare D for hopping
#should have 2.0 factor in khop Sept 14 2016
kon = koff/R  # (per second) 0.0 gives SD (singular limit)  (0.001 gives 10/s)  R=koff/kon

# density 1 at 0
# density 0 at width (which isn't a site)

print("R:",R," D0:",D0," dx:",dx," kon:",kon," koff:",koff," khop:", khop)
print("tmax:",tmax," nrun:",nrun," name:",name," seed:",seed," width:",width," arate:",arate)

#
mt = np.arange(0,width)
mtScaled = mt/(width-1)
x = np.zeros(width,int)  # 0 if empty, 1 if particle

asite = np.zeros(width,int) # 1 if acetylated, 0 if not
Nbound=Nfree=Nacetyl= 0	# number of particles
bound = np.zeros_like(x) # 1 if bound, 0 if not

x[0] = 1 #Initial Condition
rightOut = 0
leftIn = 0

def get_nth_index(lst,item,n):
    c = count()
    return next(i for i, j in enumerate(x) if j==item and next(c) == n)

def hop(pos):

	if ran.random()<0.5:        # If less than .5 than particles moves left
		if pos>0 and x[pos-1]==0: 	# pos>0 means it can move to the left, pos-1 being 0 means nothing to its left so it can move left
			x[pos] = 0
			x[pos-1]=1
			free[this]=pos-1
	else:
		
        if pos==width-1:        # If we're at the right most boundary the particles leaves out the right
            rightOut += 1		# leave to right
			x[pos]= 0
		
        elif x[pos+1]==0:	# not blocking to right
			x[pos] = 0
			x[pos+1] = 1
    return x
def acetylate(pos):
    if ran.random()<arate*dt and asite[pos]<am:
        asite[pos]+=1
        Nacetyl+=1

def progress(pos):
    if bound[pos]==0 and ran.random()*totrate<kon: #If unbound, and meets binding probability, become bound
        bound[pos]=1
        
    elif bound[pos]==0 and ran.random()*totrate<(kon+khop): #Else if unbound, and meets hopping probability
        x = hop(pos)                                        #Updates the x array because hopping occured
        
    elif bound[pos]==1 and ran.random()*totrate<koff: #Else if bound, and meets unbinding probability, unbind the particle
        bound[pos]=0
    return x,bound
TFACTOR = 2.0
### START EVOLUTION 
telapsed = 0
scale = 0
nexttime = tmax #1.0  # for next average

acetylOccupancy=np.zeros(width+1)
asite = np.zeros(width,int)
	
while telapsed<tmax:


    if x[0] == 0:  	# add to left at infinite rate
        x[0]=1
        leftIn+=1
        bound[0] = 0
    Nfree = sum(bound==0)
    Nbound = sum(bound==1)
    
    totrate = Nfree*(kon+khop)+Nbound*koff #
    dt  = -1.0/totrate*m.log(1.0-ran.random())   # kmc dt
    arate = 100*dt/totrate		 
    
    acetylOccupancy[0:width] += asite*dt
    
    telapsed += dt
    pos = get_nth_index(x,1,ran.randrange(sum(x)))
    x,bound = progress(pos)
