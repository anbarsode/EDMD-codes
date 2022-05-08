import numpy as np
import matplotlib.pyplot as plt

Nsteps = 11
fname_base = 'epd_LJ_Nsteps_%d' % Nsteps
f = open(fname_base + ".txt", "r")
while f.readline()[:3] != 'u b': True;
data = np.loadtxt(f)
f.close()

plt.subplot(131)
plt.imshow(np.reshape(data[:,7],(100,100)), origin='lower', extent=[0,1,0,10], aspect='auto')
plt.xlabel('b')
plt.ylabel('u')
plt.title(r'$\Delta T (\%)$')
plt.clim([-100,100])

plt.subplot(132)
plt.imshow(np.reshape(data[:,10],(100,100)), origin='lower', extent=[0,1,0,10], aspect='auto')
plt.title(r'$\Delta \theta (\%)$')
plt.xticks([])
plt.yticks([])
plt.clim([-100,100])

plt.subplot(133)
plt.imshow(np.reshape(data[:,4],(100,100)), origin='lower', extent=[0,1,0,10], aspect='auto')
plt.title(r'$\Delta r_{min} (\%)$')
plt.xticks([])
plt.yticks([])
plt.clim([-100,100])
plt.colorbar()

plt.suptitle('LJcut3.0, Nsteps=%d' % Nsteps)
plt.subplots_adjust(wspace=0.01, left=0.08)
plt.savefig(fname_base + ".png")
plt.show()
