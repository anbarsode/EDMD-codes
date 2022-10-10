#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <fstream>
#include <chrono>

#include "RandomICGen.cpp"
#include "Potential.cpp"
#include "Trajectory.cpp"
#include "Optimizer.cpp"

int main(int argc, char **argv)
{
    //std::ofstream f("../Data/random_nos.txt", std::ios::out);
    
    std::vector<double> params;
    params.push_back(0);
    params.push_back(3);
    params.push_back(0.5);
    
    RNG<double> uDist(3, params);
    RNG<double> bDist(4, params);
    
    std::random_device rd;
    std::mt19937_64 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    //std::mt19937_64 rng(rd());
    
    for(int Nsteps = 40; Nsteps < 100; Nsteps+=10)
    {
        Potential<double> V;
        double m=1.0, dr = 1e-6;
        V.Nsteps = Nsteps;
        for(int i=0; i<V.Nsteps; i++) V.Rsteps.push_back(V.Rcut - i*V.Rcut/(V.Nsteps-1));
        V.Rsteps[V.Nsteps-1] += dr;
        V.Vsteps.push_back(V.Vcut);
        //for(int i=1; i<V.Nsteps; i++) V.Vsteps.push_back(V.Continuous(V.Rsteps[i])); //left hand side
        for(int i=1; i<V.Nsteps; i++) V.Vsteps.push_back(V.Continuous(V.Rsteps[i-1])); //right hand side
        //for(int i=1; i<V.Nsteps; i++) V.Vsteps.push_back(V.Continuous(0.5*(V.Rsteps[i]+V.Rsteps[i-1]))); //midpoint
        
        Optimizer<double> opt(V, 3, params, 4, params, m, dr);
        //std::cout << Nsteps << " " << opt.ComputeLoss1(rng, 1e4) << std::endl; // behaves somewhat as expected
        std::cout << Nsteps << " " << opt.ComputeLoss2(rng,1e4) << std::endl; // behaves somewhat as expected
        //std::cout << Nsteps << " " << opt.ComputeLoss3(rng) << std::endl; // doesn't vary with V.Nsteps
        //std::cout << Nsteps << " " << opt.ComputeLoss4(rng) << std::endl; // wild fluctuations vs V.Nsteps
    }
    
    //f.close();
}
