# EDMD-codes
WIP. Codes for improving event driven simulations.  
  
Finding better ways to discretize continuous potentials in order to reduce computational cost in (event-driven) molecular dynamics simulations  
  
Obviously using more number of discrete steps gives trajectories that are closer to the numerically "exact" trajectories. The job is to get accurate trajectories with as few steps as possible. The discrete trajectories in the following figures were generated using a simple discretizaton scheme:  
* Uniform radial spacing between 0 to Rcut for step locations
* Value of continuous potential at $R_1$ for step potential between $R_1$ and $R_2$(>$R_1$)
<img src="https://github.com/anbarsode/EDMD-codes/blob/017ac107c3133ee54c9dab0532e99ec5dcaac773/LJcut3.0_Nsteps_101.png" alt="LJcut3.0_Nsteps_101" width="320"/>
  
Here are some example trajectories for the Lennard-Jones potential (cut off at r=3.0)  
(colors indicate how much discrete trajectory lags behind or is ahead of the continuous trajectory in terms of percentage error)  
<img src="https://github.com/anbarsode/EDMD-codes/blob/9596d1c4b43cf8a2ffed69eb8c3785c39ded736d/traj_Nsteps_101.png" alt="traj_Nsteps_101" width="320"/> <img src="https://github.com/anbarsode/EDMD-codes/blob/9596d1c4b43cf8a2ffed69eb8c3785c39ded736d/traj_Nsteps_11.png" alt="traj_Nsteps_11" width="320"/>  

The following figures show the percentage errors in discrete trajectories' times and angles with respect to numerically computed continuous trajectories for various impact parameters and speeds.  
$T$ is the time taken by the particle to escape the potential (reach Rcut again)  
$\theta$ is the angle which it moves during the same  
$r_{min}$ is the distance of closest approach  
$\Delta$ denotes the error between discrete and continuous trajectories  
$b$ is the impact parameter  
$u$ is the speed at infinity (Rcut)  
<img src="https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_1001.png" alt="traj_Nsteps_1001" width="320"/> <img src="https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_101.png" alt="traj_Nsteps_101" width="320"/> <img src="https://github.com/anbarsode/EDMD-codes/blob/6c910732b88c10441df52ac87b5e6b1c5a445d2f/epd_LJ_Nsteps_11.png" alt="traj_Nsteps_11" width="320"/>  
  
Notice how worse the percentage errors get when one reduces the number of discountinuos steps  
