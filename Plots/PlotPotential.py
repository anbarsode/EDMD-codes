import numpy as np
import matplotlib.pyplot as plt

Rsteps = '3 2.7 2.4 2.1 1.8 1.5 1.2 0.9 0.6 0.3 1e-006 '
Vsteps = '0 -0.0018626 -0.0036239 -0.00764794 -0.0178575 -0.0471002 -0.137906 -0.189378 25.9506 14382.3 7.70695e+009'
Rsteps = np.array([float(a) for a in Rsteps.split()])
Vsteps = np.array([float(a) for a in Vsteps.split()])
Nsteps = Rsteps.shape[0]

r = np.arange(1e-6,3.0,0.001)
V = 1/r**12.0 - 1/r**6.0

Vdisc = V.copy()
for i in range(len(r)):
    idx = np.argmin(Rsteps>r[i])
    Vdisc[i] = Vsteps[idx]
    
plt.plot(r,V,'k',label='Continuous')
plt.plot(r,Vdisc,'r',label='Discrete')
plt.ylim([1.1*V.min(),-2*V.min()])
plt.xlim([r[np.argmin(V>-10*V.min())]-0.1,0.1+r.max()])
plt.xlabel('r')
plt.ylabel('V(r)')
plt.title('LJcut3.0, Nsteps=%d'%Nsteps)
plt.legend()
plt.savefig('LJcut3.0_Nsteps_%d.png'%Nsteps)
plt.show()
