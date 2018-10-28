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

#SliceTimes = SliceTimes[3:] # This line is for that one broken simulation that I did, oct 19 2018

density = np.zeros_like(np.load('DENS_0.npy'))
acetylation = np.zeros_like(density)
NTot = np.zeros(np.size(TotTimes))
ATot = np.zeros_like(NTot)

densityArray = []
acetylationArray = []
NTotArray = []
ATotArray = []

flux = np.zeros(2)

for directory in directories:
        if '.npy' in directory:
                temp = np.load(directory)

        if 'DENS' in directory:
                density += temp
                densityArray.append(temp)

        elif 'ACETYL' in directory:
                acetylation += temp
                acetylationArray.append(temp)

        elif 'NTOT' in directory: 
                NTot += temp
                NTotArray.append(temp)

        elif 'ATOT' in directory:
                ATot += temp
                ATotArray.append(temp)
        elif 'FLUX' in directory:
                flux += temp
density2 = density / tubuleSims
acetylation2 = acetylation / tubuleSims
NTot /= tubuleSims
ATot /= tubuleSims
flux /= tubuleSims

os.chdir(absPath)

print('Average Inward Flux: %d , Average Outward Flux: %d'%(flux[0],flux[1]))

# Min and Max Error Bars Loops
minNTot = NTotArray[0]
maxNTot = NTotArray[0]
minATot = ATotArray[0]
maxATot = ATotArray[0]

for i in range(int(tubuleSims)):
        minNTot = np.minimum(minNTot,NTotArray[i])
        maxNTot = np.maximum(maxNTot,NTotArray[i])
        minATot = np.minimum(minATot,ATotArray[i])
        maxATot = np.maximum(maxATot,ATotArray[i])

for i, value in enumerate(SliceTimes):
        plt.figure(1)
        plt.plot(mt,density2[i],label=('t = %d s'%(value)))

        plt.figure(2)
        plt.plot(mt,acetylation2[i]/am,label=('t =%d s'%(value)))

plt.figure(1)
plt.xlabel('x/L')
plt.ylabel('Average Occupation')
plt.legend()

plt.figure(2)
plt.xlabel('x/L')
plt.ylabel('Average Acetylation')
plt.legend()

plt.figure(3)
ax = plt.gca()
plt.plot(TotTimes,NTot)
plt.fill_between(TotTimes,minNTot,maxNTot,alpha=.3)
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('t (s)')
plt.ylabel('Average Total Enzymes in Tubule')

plt.figure(4)
ax = plt.gca()
plt.plot(TotTimes,ATot)
plt.fill_between(TotTimes,minATot,maxATot,alpha=.3)
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('t (s)')
plt.ylabel('Average Total Acetylated Sites')

plt.show()
