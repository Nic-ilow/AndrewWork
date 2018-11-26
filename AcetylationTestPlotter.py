import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import datetime
import scipy.special as scis

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
SliceTimes = np.load('SliceTimes.npy')


density = np.zeros(int(width))
acetylation = np.zeros_like(SliceTimes)

densityArray = []
acetylationArray = []

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
acetylation2 = acetylation / tubuleSims

os.chdir(absPath)

z = x/np.sqrt(4*D0*SliceTimes)

acetylationfit = (1 - np.exp(-p0*arate*SliceTimes*( ((1+2*z*z)*scis.erfc(z)) - ((2*z*np.exp(-(z*z)))/np.sqrt(np.pi)) ) ) )

plt.scatter(SliceTimes,acetylation2/am,label='Stochastic')
plt.plot(SliceTimes,acetylationfit,label='Analytic',c='r')
plt.xlabel('time (s)')
plt.ylabel('Fraction of Acetylated Sites')
plt.legend()
plt.show()
