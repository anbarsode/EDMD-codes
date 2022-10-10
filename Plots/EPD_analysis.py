import numpy as np
import matplotlib.pyplot as plt
from myplotsettings import *

Nsteps = 11
fname_base = 'epd_LJ_Nsteps_%d' % Nsteps
f = open(fname_base + ".txt", "r")
while f.readline()[:3] != 'u b': True;
data = np.loadtxt(f)
f.close()

clim = [-10,10]

plt.subplot(131)
plt.imshow(np.reshape(data[:,7],(200,100)), origin='lower', extent=[0,1,0,10], aspect='auto')
plt.xlabel('b')
plt.ylabel('u')
plt.title(r'$\Delta T (\%)$')
plt.clim(clim)

plt.subplot(132)
plt.imshow(np.reshape(data[:,10],(200,100)), origin='lower', extent=[0,1,0,10], aspect='auto')
plt.title(r'$\Delta \theta (\%)$')
plt.xticks([])
plt.yticks([])
plt.clim(clim)

plt.subplot(133)
plt.imshow(np.reshape(data[:,4],(200,100)), origin='lower', extent=[0,1,0,10], aspect='auto')
plt.title(r'$\Delta r_{min} (\%)$')
plt.xticks([])
plt.yticks([])
plt.clim(clim)
plt.colorbar()

plt.suptitle('LJcut3.0, Nsteps=%d' % Nsteps)
plt.tight_layout()
plt.subplots_adjust(wspace=0.01, left=0.1, right=0.92)
plt.savefig(fname_base + ".png")
plt.show()
