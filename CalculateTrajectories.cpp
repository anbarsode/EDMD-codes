//$ g++ CalculateTrajectories.cpp
//$ ./a.out

#include <iostream>
#include <math.h>
#include <functional>
#include <vector>
#include <fstream>
#include <chrono>

#include "Potential.cpp"
#include "Trajectory.cpp"

int main()
{
    std::ofstream cont("./cont_traj.txt", std::ios::out);
    std::ofstream disc("./disc_traj.txt", std::ios::out);
    
    Potential<double> V;
    double u=0.5, b=0.2, m=1.0, dr = 1e-6;
    b *= V.Rcut;
    
    V.Nsteps = 101;
    for(int i=0; i<V.Nsteps; i++) V.Rsteps.push_back(V.Rcut - i*V.Rcut/(V.Nsteps-1));
    V.Rsteps[V.Nsteps-1] += dr;
    for(int i=0; i<V.Nsteps; i++) V.Vsteps.push_back(V.Continuous(V.Rsteps[i]));
    
    //for(int i=0; i<V.Nsteps; i++) std::cout << V.Rsteps[i] << " " << V.Vsteps[i] << std::endl;
    
    
    auto start = std::chrono::high_resolution_clock::now();
    
    Trajectory<double> traj(V, u, b, m, dr);
    traj.ComputeTrajectories();
    
    for(int i=0; i<traj.ContR.size(); i++) cont << std::defaultfloat << traj.ContR[i] << " " << std::scientific << traj.ContDt[i] << " " << traj.ContDth[i] << std::endl;
    for(int i=0; i<traj.DiscR.size(); i++) disc << std::defaultfloat << traj.DiscR[i] << " " << std::scientific << traj.DiscDt[i] << " " << traj.DiscDth[i] << std::endl;
    
    cont.close();
    disc.close();
    
    auto elap = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "Elapsed time (s) : " << elap.count()/1e6 << std::endl;
    return 0;
}
