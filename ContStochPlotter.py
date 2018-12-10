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
parser.add_argument('-sd','--sd',default=now.day)
parser.add_argument('-sm','--sm',default=now.month)
parser.add_argument('-sy','--sy',default=now.year)
parser.add_argument('-shour','--shour',default=now.hour)
parser.add_argument('-smin','--smin',default=now.minute)
parser.add_argument('-cd','--cd',default=now.day)
parser.add_argument('-cm','--cm',default=now.month)
parser.add_argument('-cy','--cy',default=now.year)
parser.add_argument('-chour','--chour',default=now.hour)
parser.add_argument('-cmin','--cmin',default=now.minute)
args = parser.parse_args()

sday = int(args.sd)
smonth = int(args.sm)
syear = int(args.sy)
shour = int(args.shour)
sminute = int(args.smin)

cday = int(args.cd)
cmonth = int(args.cm)
cyear = int(args.cy)
chour = int(args.chour)
cminute = int(args.cmin)

stochDataPath = 'Sim_{0}_{1}_{2}_{3}_{4}'.format(sday,smonth,syear,shour,sminute)
ContDataPath = 'Sim_{0}_{1}_{2}_{3}_{4}'.format(cday,cmonth,cyear,chour,cminute)

absPath = os.getcwd()

os.chdir('StochasticData')
os.chdir(stochDataPath)

density = np.zeros_like(np.load('DENS_0.npy'))
StochDirectories = os.listdir(os.getcwd())
StochSimInfo = np.loadtxt('SIMINFO')
TotTimes = np.load('TotTimes.npy')
SliceTimes = np.load('SliceTimes.npy')
width,tmax,tubuleSims,arate,am,koff,kon = StochSimInfo


os.chdir(absPath)
os.chdir('ContinuousData')
os.chdir(ContDataPath)

ContDirectories = os.listdir(os.getcwd())
ContSimInfo = np.loadtxt('siminfo')
dx,width,kon,koff,dt,arate,tmax = ContSimInfo



mt=np.arange(0,width)/width
x = mt*width*dx
D0 = 2.7e5
p0 = 1.0

acetylation = np.zeros_like(density)
NTot = np.zeros(np.size(TotTimes))
ATot = np.zeros_like(NTot)

ContDensity = np.zeros_like(np.load('DENS.npy'))
ContAcetylation = np.zeros_like(ContDensity)
ContNTot = np.zeros_like(NTot)
ContATot = np.zeros_like(ATot)

densityArray = []
acetylationArray = []
NTotArray = []
ATotArray = []

flux = np.zeros(2)

for directory in ContDirectories: ##Loading in Continuous Data
        if '.npy' in directory:
                temp = np.load(directory)

        if 'DENS' in directory:
                ContDensity += temp

        elif 'ACETYL' in directory:
                ContAcetylation += temp

        elif 'NTOT' in directory:
                ContNTot += temp

        elif 'ATOT' in directory:
                ContATot += temp
os.chdir(absPath)
os.chdir('StochasticData')
os.chdir(stochDataPath)
for directory in StochDirectories: ##Loading in Stochastic Data
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

def aAnalytic(t):
        global p0,arate,D0,x
        z = x/np.sqrt(4*D0*t)
        afit = 1 - np.exp(-p0*arate*t*( ((1+2*z*z)*scis.erfc(z)) - ((2*z*np.exp(-(z*z)))/np.sqrt(np.pi)) ) )
        return afit

def pAnalytic(t):
        z = x/np.sqrt(4*D0*t)
        pfit = scis.erfc(z) * p0
        return pfit

for i, value in enumerate(SliceTimes):
        
        plt.figure(1)
        plt.plot(mt,density2[i],label=('Stochastic Model t = %.2f s'%(value)))
        plt.plot(mt,ContDensity[i],ls='dashed',label=('Continuous Model t= %.2f s'%(value)))
        plt.figure(2)
        plt.plot(mt,acetylation2[i]/am,label=('Stochastic Model t = %.2f s'%(value)))
        plt.plot(mt,ContAcetylation[i],ls='dashed',label=('Continuous Model t = %.2f s'%(value)))

plt.figure(3)

NtotAnalytic = p0*np.sqrt(4*D0/(dx*dx)*TotTimes/np.pi)


plt.figure(1)
plt.xlabel('x/L')
plt.ylabel('Density')
plt.legend()

plt.figure(2)
plt.xlabel('x/L')
plt.ylabel('Acetylation')
plt.legend()

plt.figure(3)
ax = plt.gca()
plt.plot(TotTimes,NTot,label='Stochastic Model')
plt.plot(TotTimes,ContNTot,ls='dashed',label='Continuous Model')
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('t (s)')
plt.ylabel('Total Enzymes in Tubule')
plt.legend()

ATotFit = []
for t in TotTimes:
        ATotFit.append(sum(aAnalytic(t)))

plt.figure(4)
ax = plt.gca()
plt.plot(TotTimes,ATot,label='Stochastic Model')
plt.plot(TotTimes,ContATot,ls='dashed',label='Continuous Model')
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('t (s)')
plt.ylabel('Total Acetylated Sites')
plt.legend()

plt.show()
