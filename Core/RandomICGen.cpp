template <class T = double>
class RNG
{
    public:
        int PDF;
        std::vector<T> params;
        T RootKTBymu, MinVal, MaxVal;
        T x, y, z;
        const T PI = atan(1.0) * 4.0;
        
        RNG(int PDF_, std::vector<T> params_)
        {
            PDF = PDF_;
            params = params_;
        }
        
        T Dist(std::mt19937_64& rng)
        {
            switch(PDF)
            {
                case 0: // Continuous uniform distribution
                {
                    if(params[0] >= params[1]) std::cout << "Warning: lower_limit >= upper_limit in uniform distribution" << std::endl;
                    std::uniform_real_distribution<T> cu(params[0], params[1]);
                    return cu(rng);
                    break;
                }
                
                case 1: // 1D relative Maxwell-Boltzmann
                {
                    //if(params[2] <= 0) std::cout << "Warning: Temperature <= 0 in Maxwell-Boltzmann distribution" << std::endl;
                    RootKTBymu = params[2];
                    
                    if(params[0] < 0) MinVal = 0;
                    else MinVal = params[0];
                    
                    if(params[1] < 0) MaxVal = 6 * RootKTBymu; // roughly 1 in a million chance of missing higher values
                    else MaxVal = params[1];
                    
                    T RootTwoKTByPimu = sqrt(2.0 / PI) / RootKTBymu;
                    std::uniform_real_distribution<T> cu(MinVal, MaxVal);
                    std::uniform_real_distribution<T> cu_base(0, 1);
                    do
                    {
                        x = cu(rng);
                        y = cu_base(rng);
                        z = RootTwoKTByPimu * exp(-0.5 * pow(x / RootKTBymu, 2.0));
                    }while(y > z);
                    return x;
                    break;
                }
                
                case 2: // 2D relative Maxwell-Boltzmann
                {
                    //if(params[2] <= 0) std::cout << "Warning: Temperature <= 0 in Maxwell-Boltzmann distribution" << std::endl;
                    RootKTBymu = params[2];
                    
                    if(params[0] < 0) MinVal = 0;
                    else MinVal = params[0];
                    
                    if(params[1] < 0) MaxVal = 6 * RootKTBymu; // roughly 1 in a million chance of missing higher values
                    else MaxVal = params[1];
                    
                    T muByKT = 1 / RootKTBymu / RootKTBymu;
                    std::uniform_real_distribution<T> cu(MinVal, MaxVal);
                    std::uniform_real_distribution<T> cu_base(0, 1);
                    do
                    {
                        x = cu(rng);
                        y = cu_base(rng);
                        z = muByKT * x * exp(-0.5 * pow(x / RootKTBymu, 2.0));
                    }while(y > z);
                    return x;
                    break;
                }
                
                case 3: // 3D relative Maxwell-Boltzmann
                {
                    //if(params[2] <= 0) std::cout << "Warning: Temperature <= 0 in Maxwell-Boltzmann distribution" << std::endl;
                    RootKTBymu = params[2];
                    
                    if(params[0] < 0) MinVal = 0;
                    else MinVal = params[0];
                    
                    if(params[1] < 0) MaxVal = 6 * RootKTBymu; // roughly 1 in a million chance of missing higher values
                    else MaxVal = params[1];
                    
                    T RootTwoByPiByaCube = sqrt(2.0 / PI) / RootKTBymu / RootKTBymu / RootKTBymu;
                    std::uniform_real_distribution<T> cu(MinVal, MaxVal);
                    std::uniform_real_distribution<T> cu_base(0, 1);
                    do
                    {
                        x = cu(rng);
                        y = cu_base(rng);
                        z = RootTwoByPiByaCube * x * x * exp(-0.5 * pow(x / RootKTBymu, 2.0));
                    }while(y > z);
                    return x;
                    break;
                }
                
                case 4: // Continuous linear distribution from MinVal to MaxVal
                {
                    //if(params[0] >= params[1]) std::cout << "Warning: lower_limit >= upper_limit in linear distribution" << std::endl;
                    MinVal = params[0];
                    MaxVal = params[1];
                    std::uniform_real_distribution<T> cu(MinVal, MaxVal);
                    std::uniform_real_distribution<T> cu_base(0, 2 / (MaxVal - MinVal));
                    do
                    {
                        x = cu(rng);
                        y = cu_base(rng);
                        z = 2 * (x - MinVal) / (MaxVal - MinVal) / (MaxVal - MinVal);
                    }while(y > z);
                    return x;
                    break;
                }
            }
        }
};

/*
template <class T = double>
class ICGen
{
    public:
        std::string fname;
        RNG uDist, bDist;
        
        ICGen(std::string fname_)
        {
            fname = fname_;
        }
        
        void ReadDistributions(void)
        {
            // Read distributions from fname
            ;
        }
        
        T* uRNG(int N, std::mt19937_64& rng)
};
*/