//$ g++ EndPointDeviation.cpp
//$ ./a.out

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <chrono>

#include "Potential.cpp"
#include "Trajectory.cpp"

int main()
{
    std::ofstream epd("./epd_LJ_Nsteps_11.txt", std::ios::out);
    
    double m = 1.0, dr = 1e-6;
    double umin = 1e-6, umax = 10.0, du = 0.1;
    double bmin = 0.0, bmax = 1.0, db = 0.01;
    
    Potential<double> V;    
    V.Nsteps = 11;
    for(int i=0; i<V.Nsteps; i++) V.Rsteps.push_back(V.Rcut - i*V.Rcut/(V.Nsteps-1));
    V.Rsteps[V.Nsteps-1] += dr;
    for(int i=0; i<V.Nsteps; i++) V.Vsteps.push_back(V.Continuous(V.Rsteps[i]));
    
    epd << "Metadata" << std::endl;
    epd << "Potential : LJ cut, 1.0/r12 - 1.0 / r6" << std::endl;
    epd << "Nsteps = " << V.Nsteps << std::endl;
    epd << "Rsteps :" << std::endl;
    for(int i=0; i<V.Nsteps; i++) epd << V.Rsteps[i] << " ";
    epd << std::endl;
    epd << "Vsteps :" << std::endl;
    for(int i=0; i<V.Nsteps; i++) epd << V.Vsteps[i] << " ";
    epd << std::endl;
    epd << "\nm = " << std::defaultfloat << m << std::endl;
    epd << "dr = " << std::defaultfloat << dr << std::endl;
    epd << "u = " << std::defaultfloat << umin << ":" << umax << ":" << du << std::endl;
    epd << "b = " << std::defaultfloat << bmin << ":" << bmax << ":" << db << std::endl;
    epd << "\nEnd Point Values" << std::endl;
    epd << "u b ContRmin DiscRmin %ERmin ContDtEnd DiscDtEnd %EDt ContDthEnd DiscDthEnd %EDth" << std::endl;
        
    auto start = std::chrono::high_resolution_clock::now();
    for(double u = umin; u < umax; u += du)
    {
        for(double b = bmin; b < bmax; b += db)
        {
            Trajectory<double> traj(V, u, b * V.Rcut, m, dr);
            traj.EndPointValues();
            epd << std::defaultfloat << u << " " << b << std::scientific
            << " " << traj.ContRmin << " " << traj.DiscRmin << " " << traj.DiscRmin / traj.ContRmin * 100 - 100
            << " " << traj.ContDtEnd << " " << traj.DiscDtEnd << " " << traj.DiscDtEnd / traj.ContDtEnd * 100 - 100
            << " " << traj.ContDthEnd << " " << traj.DiscDthEnd << " " << traj.DiscDthEnd / traj.ContDthEnd * 100 - 100
            << std::endl;
        }
        auto elap = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
        std::cout << (u - umin) / (umax - umin) * 100.0 << "% done. Elapsed time (s) : " << elap.count()/1e6 << std::endl;
        //std::cout << (u - umin) / (umax - umin) * 100.0 << "% done. Elapsed time (s) : " << elap.count()/1e6 << "\r";
    }
    epd.close();
    return 0;
}
