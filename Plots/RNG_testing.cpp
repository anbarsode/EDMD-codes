#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <fstream>

#include "RandomICGen.cpp"

int main(int argc, char **argv)
{
    std::ofstream f("./random_nos.txt", std::ios::out);
    
    std::vector<double> params;
    params.push_back(0);
    params.push_back(6);
    params.push_back(1);
    
    RNG<double> cu(0, params);
    RNG<double> mb1(1, params);
    RNG<double> mb2(2, params);
    RNG<double> mb3(3, params);
    RNG<double> cl(4, params);
    
    std::random_device rd;
    std::mt19937_64 rng(rd());
    
    for(long i=0; i<100000; i++)
        f << cu.Dist(rng) << " " << mb1.Dist(rng) << " " << mb2.Dist(rng) << " " << mb3.Dist(rng) << " " << cl.Dist(rng) << std::endl;
    
    f.close();
}