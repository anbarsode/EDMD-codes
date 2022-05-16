import numpy as np
import matplotlib.pyplot as plt

a = 1.0
x = np.linspace(0,6,1000)
cu = np.ones_like(x) / (x.max() - x.min())
cl = x * 2 / (x.max() - x.min())**2.0
m1 = (2.0/np.pi)**0.5 / a*np.exp(-x**2.0/2.0/a**2.0)
m2 = x / a**2.0 * np.exp(-x**2.0/2.0/a**2.0)
m3 = (2.0/np.pi)**0.5 * x**2.0 / a**3.0 * np.exp(-x**2.0/2.0/a**2.0)

rn = np.loadtxt('./random_nos.txt')
nbins = 50

plt.subplot(221)
plt.hist(rn[:,0], nbins, color='r', density=True, alpha=0.5)
plt.hist(rn[:,-1], nbins, color='b', density=True, alpha=0.5)
plt.plot(x,cu,'r',label='uniform')
plt.plot(x,cl,'b',label='linear')
plt.legend()
#plt.ylim([1e-6,1])
#plt.yscale('log')

plt.subplot(222)
plt.hist(rn[:,1], nbins, color='b', density=True, alpha=0.5)
plt.plot(x,m1,'b',label='MB1D')
plt.legend()
#plt.ylim([1e-6,1])
#plt.yscale('log')

plt.subplot(223)
plt.hist(rn[:,2], nbins, color='b', density=True, alpha=0.5)
plt.plot(x,m2,'b',label='MB2D')
plt.legend()
#plt.ylim([1e-6,1])
#plt.yscale('log')

plt.subplot(224)
plt.hist(rn[:,3], nbins, color='b', density=True, alpha=0.5)
plt.plot(x,m3,'b',label='MB3D')
plt.legend()
#plt.ylim([1e-6,1])
#plt.yscale('log')

plt.savefig('./RNG_testing.png')
plt.show()