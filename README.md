# EDMD-codes
WIP. Codes for improving event driven simulations.  
  
Finding better ways to discretize continuous potentials in order to reduce computational cost in (event-driven) molecular dynamics simulations  
  
Obviously using more number of discrete steps gives trajectories that are closer to the numerically "exact" trajectories. The job is to get accurate trajectories with as few steps as possible. The discrete trajectories in the following figures were generated using a simple discretizaton scheme:  
* Uniform radial spacing between 0 to Rcut for step locations
* Value of continuous potential at $(R_1+R_2)/2$ for step potential between $R_1$ and $R_2(>R_1)$  
 
<img src="https://github.com/anbarsode/EDMD-codes/blob/8d9cec1d31de3a8db353b5f8033e6279ded1be23/Plots/LJcut3.0_Nsteps_11.png" alt="LJcut3.0_Nsteps_11" width="320"/>
  
Here are some example trajectories for the Lennard-Jones potential (cut off at r=3.0)  
(colors indicate how much discrete trajectory lags behind or is ahead of the continuous trajectory in terms of percentage error)  
<img src="https://github.com/anbarsode/EDMD-codes/blob/8d9cec1d31de3a8db353b5f8033e6279ded1be23/Plots/traj_Nsteps_11.png" alt="traj_Nsteps_11" width="320"/> <img src="https://github.com/anbarsode/EDMD-codes/blob/8d9cec1d31de3a8db353b5f8033e6279ded1be23/Plots/traj_Nsteps_101.png" alt="traj_Nsteps_101" width="320"/>  

The following figures show the percentage errors in discrete trajectories' times and angles with respect to numerically computed continuous trajectories for various impact parameters and speeds.  
$T$ is the time taken by the particle to escape the potential (reach Rcut again)  
$\theta$ is the angle which it moves during the same  
$r_{min}$ is the distance of closest approach  
$\Delta$ denotes the error between discrete and continuous trajectories  
$b$ is the impact parameter  
$u$ is the speed at infinity (Rcut)  
<img src="https://github.com/anbarsode/EDMD-codes/blob/008d2f118f94e59b74638bfa90a0ba183d654ab4/Plots/epd_LJ_Nsteps_11.png" alt="epd_LJ_Nsteps_11" width="320"/> <img src="https://github.com/anbarsode/EDMD-codes/blob/008d2f118f94e59b74638bfa90a0ba183d654ab4/Plots/epd_LJ_Nsteps_101.png" alt="epd_LJ_Nsteps_101" width="320"/>  
  
Tests on random number generation based on a given distribution (to be used later during optimization)  
<img src="https://github.com/anbarsode/EDMD-codes/blob/4eca4377eade6885b82e2808f158e3cd8a3dc27e/Plots/RNG_testing.png" alt="RNG_testing" width="640"/>
