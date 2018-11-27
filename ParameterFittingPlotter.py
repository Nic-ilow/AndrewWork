import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import datetime
import scipy.special as scis
import scipy.optimize as scio

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

dx = 7.0
x = mt*width*dx
D0 = 2.7e5
p0 = 1.0
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
fudge = am/(dx+1.0)

def aAnalytic(x,Beta):
        global p0,arate,D0,value
        z = x/np.sqrt(4*D0*value)
        afit = 1 - np.exp(-p0*Beta*arate*value*( ((1+2*z*z)*scis.erfc(z)) - ((2*z*np.exp(-(z*z)))/np.sqrt(np.pi)) ) )
        return afit

def pAnalytic(t):
        z = x/np.sqrt(4*D0*t)
        pfit = scis.erfc(z) * p0
        return pfit
fitParam = []
for i, value in enumerate(SliceTimes):
        if kon==0:
                temp = scio.curve_fit(aAnalytic,x,acetylation2[i]/am)
                fitParam.append(temp[0])
                acetylationfit = aAnalytic(x,1)
                pfit = pAnalytic(value)
                plt.figure(1)
                plt.plot(mt,pfit,ls='dashed',lw=4,label=('Analytic t=%.2f s'%(value)))
                plt.figure(2)
                plt.plot(mt,acetylationfit,ls='dashed',lw=4,label=('Analytic t=%.2f s'%(value)))
                plt.plot(mt,aAnalytic(x,temp[0]),lw=3,label=('Optimized Fudge Factor =%.3f'%(temp[0])))
        plt.figure(1)
        plt.plot(mt,density2[i],label=('t = %.2f s'%(value)))
        plt.figure(2)
        plt.plot(mt,acetylation2[i]/am,label=('t =%.2f s'%(value)))
if kon==0:
        ax = plt.gca()
        plt.figure(3) 
        NtotAnalytic = p0*np.sqrt(4*D0/(dx*dx)*TotTimes/np.pi)
        plt.plot(TotTimes,NtotAnalytic,c='r',ls='dashed',lw=4,label='Analytic')

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
plt.plot(TotTimes,NTot,label='Stochastic')
plt.fill_between(TotTimes,minNTot,maxNTot,alpha=.3)
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('t (s)')
plt.ylabel('Average Total Enzymes in Tubule')
plt.legend()

ATotFit = []
for value in TotTimes:
        ATotFit.append(sum(aAnalytic(x,1)))

plt.figure(4)
ax = plt.gca()
plt.plot(TotTimes,ATot,label='Stochastic')
plt.plot(TotTimes,ATotFit,label='Analytic')
plt.fill_between(TotTimes,minATot,maxATot,alpha=.3)
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('t (s)')
plt.ylabel('Average Total Acetylated Sites')
plt.legend()

plt.show()
