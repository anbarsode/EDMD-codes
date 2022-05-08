import numpy as np
import matplotlib.pyplot as plt

Nsteps = 101
Rcut = 3.0
u = 2.0
b = 0.3
m = 1.0
theta0 = np.arcsin(b)

contR, contDt, contDth = np.loadtxt('./cont_traj.txt',unpack=True)
contDt = contDt.cumsum() # remove or add [:-1] after contDt if you get an error
contDth = contDth.cumsum()

discR, discDt, discDth = np.loadtxt('./disc_traj.txt',unpack=True)
discDt = discDt.cumsum()
discDth = discDth.cumsum()

# scipy traj
R = np.linspace(Rcut,contR.min()-0.1,10000)
def pot(r):
    return 1.0/r**12.0 - 1.0/r**6.0

def ydot(r,y,u,b,m):
    dtdr = -1/np.sqrt(u**2.0*(1-b**2.0/r**2.0)-2.0/m*pot(r))
    return [dtdr, u*b*dtdr/r/r]
    
from scipy.integrate import solve_ivp
sol = solve_ivp(ydot, [R[0], R[-1]], [0, theta0], args=(u,b*Rcut,m), t_eval=R)
xs = sol.t * np.cos(sol.y[1])
ys = sol.t * np.sin(sol.y[1])

print('Rmin: scipy %f, cont %f, disc %f' % (sol.t.min(), contR.min(), discR.min()))
plt.plot(sol.t,sol.y[1])
plt.plot(contR,contDth + theta0)
plt.plot(discR,discDth + theta0)
plt.show()

xc = contR * np.cos(contDth + theta0)
yc = contR * np.sin(contDth + theta0)
xd = discR * np.cos(discDth + theta0)
yd = discR * np.sin(discDth + theta0)

plt.plot(xs,ys,'r',alpha=0.5,label='Scipy')
plt.plot(xc,yc,'*k',alpha=0.7,ms=3,label='Continuous')
plt.scatter(xd,yd,c=discDt/contDt*100-100,s=20,label='Discrete')

cbar = plt.colorbar()
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$\Delta T (\%)$', rotation=0)

'''
plt.axhline(y=0,linestyle='--',alpha=0.3,color='b')
plt.axvline(x=0,linestyle='--',alpha=0.3,color='b')
plt.plot(np.arange(-Rcut,Rcut+0.01,0.01),np.sqrt(Rcut**2.0 - np.arange(-Rcut,Rcut+0.01,0.01)**2.0),linestyle='--',alpha=0.3,color='b')
plt.axis('equal')
xlim = [np.min([xs.min(),xc.min(),xd.min()]) - Rcut / 10.0, np.max([xs.max(),xc.max(),xd.max()]) + Rcut / 10.0]
ylim = [np.min([ys.min(),yc.min(),yd.min()]) - Rcut / 10.0, np.max([ys.max(),yc.max(),yd.max()]) + Rcut / 10.0]
plt.xlim(xlim)
plt.ylim(ylim)
'''

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('LJcut3.0, u = %.2f, b = %.2f, Nsteps=%d' % (u,b,Nsteps))
plt.savefig('./traj_Nsteps_%d.png' % Nsteps)
plt.show()
