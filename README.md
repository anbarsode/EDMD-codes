# EDMD-codes
WIP. Codes for improving event driven simulations.  
  
Finding better ways to discretize continuous potentials in order to reduce computational cost in (event-driven) molecular dynamics simulations  
  
Obviously using more number of discrete steps gives trajectories that are closer to the numerically "exact" trajectories. The job is to get accurate trajectories with as few steps as possible.  
For now, check out how worse the percentage errors get when one reduces the number of discountinuos steps  
$T$ is the time taken by the particle to escape the potential (reach Rcut again)  
$\theta$ is the angle which it moves during the same  
$r_{min}$ is the distance of closest approach  
$\Delta$ denotes the error between discrete and continuous trajectories  
$b$ is the impact parameter  
$u$ is the speed at infinity (Rcut)  

![epd_LJ_Nsteps_1001](https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_1001.png)
![epd_LJ_Nsteps_101](https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_101.png)
![epd_LJ_Nsteps_11](https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_11.png)
