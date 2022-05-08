import numpy as np
import matplotlib.pyplot as plt

Nsteps = 101
u = 0.5
b = 0.2
theta0 = np.arcsin(b)

contR, contDt, contDth = np.loadtxt('./cont_traj.txt',unpack=True)
contR = np.concatenate((contR, contR[-2::-1]))
contDt = np.concatenate((contDt[:-1], contDt[-2:0:-1])).cumsum()
contDth = np.concatenate((contDth, contDth[-1:0:-1])).cumsum()

discR, discDt, discDth = np.loadtxt('./disc_traj.txt',unpack=True)
discR = np.concatenate((discR, discR[-2::-1]))
discDt = np.concatenate((discDt, discDt[-1:0:-1])).cumsum()
discDth = np.concatenate((discDth, discDth[-1:0:-1])).cumsum()

xc = contR * np.cos(contDth + theta0)
yc = contR * np.sin(contDth + theta0)
xd = discR * np.cos(discDth + theta0)
yd = discR * np.sin(discDth + theta0)

plt.plot(xc,yc,'*k',alpha=0.7,ms=3,label='Continuous')
plt.scatter(xd,yd,c=discDt/contDt*100-100,s=7,label='Discrete')
cbar = plt.colorbar()
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$\Delta T (\%)$', rotation=0)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('LJcut3.0, u = %.2f, b = %.2f, Nsteps=%d' % (u,b,Nsteps))
plt.savefig('./traj_Nsteps_%d.png' % Nsteps)
plt.show()