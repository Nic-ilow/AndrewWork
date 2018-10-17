import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import datetime

now = datetime.datetime.now()

parser=argparse.ArgumentParser()
parser.add_argument('-d','--d',default=now.day)
parser.add_argument('-m','--m',default=now.month)
parser.add_argument('-y','--y',default=now.year)
parser.add_argument('-hour','--hour',default=now.hour)
parser.add_argument('-min','--min',default=now.minute)
args = parser.parse_args()

day = int(args.d)
month = int(args.m)
year = int(args.y)
hour = int(args.hour)
minute = int(args.min)

dataPath = 'Sim_{0}_{1}_{2}_{3}_{4}'.format(day,month,year,hour,minute)

absPath = os.getcwd()
os.chdir('Data')
os.chdir(dataPath)

directories = os.listdir(os.getcwd())

simInfo = np.loadtxt('SIMINFO')

width,tmax,tubuleSims,arate,am,koff,kon = simInfo

mt=np.arange(0,width)/width
TotTimes = np.load('TotTimes.npy')
SliceTimes = np.load('SliceTimes.npy')

density = np.zeros_like(np.load('DENS_0.npy'))
acetylation = np.zeros_like(density)
NTot = np.zeros(np.size(TotTimes))
ATot = np.zeros_like(NTot)
flux = np.zeros(2)

for directory in directories:
        
        if ('DENS') in directory:
                density += np.load(directory)

        elif 'ACETYL' in directory:
                acetylation += np.load(directory)

        elif 'NTOT' in directory: 
                NTot += np.load(directory)

        elif 'ATOT' in directory:
                ATot += np.load(directory)
        
        elif 'FLUX' in directory:
                flux += np.load(directory)

density /= tubuleSims
acetylation /= tubuleSims
NTot /= tubuleSims
ATot /= tubuleSims
flux /= tubuleSims

os.chdir(absPath)

print('Average Inward Flux: %d , Average Outward Flux: %d'%(flux[0],flux[1]))
for i, value in enumerate(SliceTimes):
        plt.figure(1)
        plt.plot(mt,density[i]/SliceTimes[i],label=('T = %d s'%(value)))

        plt.figure(2)
        plt.plot(mt,acetylation[i]/SliceTimes[i],label=('T =%d s'%(value)))

plt.figure(1)
plt.xlabel('X')
plt.ylabel('Average Density')
plt.legend()

plt.figure(2)
plt.xlabel('X')
plt.ylabel('Average Acetylation')
plt.legend()

plt.figure(3)
plt.loglog(TotTimes,NTot)
plt.xlabel('T (s)')
plt.ylabel('Average NTot')

plt.figure(4)
plt.loglog(TotTimes,ATot)
plt.xlabel('T (s)')
plt.ylabel('Average ATot')

plt.show()
