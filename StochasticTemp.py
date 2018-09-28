import math as m
import numpy as np
import matplotlib.pyplot as plt
import random as ran
import argparse
import time
import pickle
from itertools import count
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
parser.add_argument('-am','--am',default=1,help='acetyl multiplicity, e.g. how many acetyl groups per site')
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
global x
x = np.zeros(width,int)  # 0 if empty, 1 if particle

asite = np.zeros(width,int) # 1 if acetylated, 0 if not
Nbound=Nfree=Nacetyl= 0	# number of particles
free = []
bound = []

x[0] = 1 #Initial Condition
rightOut = 0
leftIn = 0

def hop(this):
        global leftIn,rightOut,Nfre
        pos = free[this]

        if ran.random()<0.5:
                if pos>0 and x[pos=1]==0:
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
def acetylate(this):
        temp = ran.randrange(Nfree+Nbound)
        if temp<Nfree:
                pos = free[temp]
                if asite[pos]<am:
                        asite[pos]+=1
        else:
                pos = bound[temp-Nfree]
                if asite[pos]<am:
                        asite[pos]+=1

