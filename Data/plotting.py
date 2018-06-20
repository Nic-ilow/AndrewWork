import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy.special as scis
dx = 7.0
tmax = 1.0
width=512
p0=1
D0 = 2.7E5
x=np.arange(width)*dx
x2 = x/width/dx
steadyState = np.linspace(p0,0,width)
fitcurve = scis.erfc(x/np.sqrt(4*D0*tmax))*steadyState/p0

#acetyl12 = np.load('acetyl1000000000000.0.npy')
acetyl1000000 = np.load('acetyl1000000.0.npy')
acetyl10 = np.load('acetyl10.0.npy')
acetyl1 = np.load('acetyl1.0.npy')
acetyl01 = np.load('acetyl0.1.npy')
acetyl001 = np.load('acetyl0.01.npy')
acetyl0001 = np.load('acetyl0.001.npy')

occupation12 = np.load('occupation1000000000000.0.npy')
occupation1000000 = np.load('occupation1000000.0.npy')
occupation10 = np.load('occupation10.0.npy')
occupation1 = np.load('occupation1.0.npy')
occupation01 = np.load('occupation0.1.npy')
occupation001 = np.load('occupation0.01.npy')
occupation0001 = np.load('occupation0.001.npy')

plt.title('tmax = 1s OCCUPATION')
plt.xlabel('normalized length')
plt.ylabel('fraction of tmax')
#plt.plot(x2,occupation12,label='pscale = 1E12')
plt.plot(x2,occupation1000000,label='pscale = 1000000')
plt.plot(x2,occupation10,label='pscale = 10.0')
plt.plot(x2,occupation1,label='pscale = 1.0')
plt.plot(x2,occupation01,label='pscale = 0.1')
plt.plot(x2,occupation001,label='pscale = 0.01')
plt.plot(x2,occupation0001,label='pscale = 0.001')
plt.plot(x2,fitcurve,label='No SFD effects Fit')
plt.legend()

plt.figure()
plt.title('tmax = 1s ACETYLATION')
plt.xlabel('normalized length')
plt.ylabel('fraction of tmax')
#plt.plot(x2,acetyl12,label='pscale = 1E12')
plt.plot(x2,acetyl1000000,label='pscale = 1000000')
plt.plot(x2,acetyl10,label='pscale = 10.0')
plt.plot(x2,acetyl1,label='pscale = 1.0')
plt.plot(x2,acetyl01,label='pscale = 0.1')
plt.plot(x2,acetyl001,label='pscale = 0.01')
plt.plot(x2,acetyl0001,label='pscale = 0.001')
plt.legend()
plt.show()

