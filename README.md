# EDMD-codes
WIP. Codes for improving event driven simulations.  
  
Finding better ways to discretize continuous potentials in order to reduce computational cost in (event-driven) molecular dynamics simulations  
  
Obviously using more number of discrete steps gives trajectories that are closer to the numerically "exact" trajectories. The job is to get accurate trajectories with as few steps as possible. The discrete trajectories in the following figures were generated using a simple discretizaton scheme:  
* Uniform radial spacing between 0 to Rcut for step locations
* Value of continuous potential at step location for step potential
  
Here are some sample trajectories for the Lennard-Jones potential (cut off at r=3.0)  
![traj_Nsteps_101](https://github.com/anbarsode/EDMD-codes/blob/9596d1c4b43cf8a2ffed69eb8c3785c39ded736d/traj_Nsteps_101.png)
![traj_Nsteps_11](https://github.com/anbarsode/EDMD-codes/blob/9596d1c4b43cf8a2ffed69eb8c3785c39ded736d/traj_Nsteps_11.png)

The following figures show the percentage errors in discrete trajectories' times and angles with respect to numerically computed continuous trajectories.  
$T$ is the time taken by the particle to escape the potential (reach Rcut again)  
$\theta$ is the angle which it moves during the same  
$r_{min}$ is the distance of closest approach  
$\Delta$ denotes the error between discrete and continuous trajectories  
$b$ is the impact parameter  
$u$ is the speed at infinity (Rcut)  
![epd_LJ_Nsteps_1001](https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_1001.png)
![epd_LJ_Nsteps_101](https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_101.png)
![epd_LJ_Nsteps_11](https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_11.png)
  
Notice how worse the percentage errors get when one reduces the number of discountinuos steps  
