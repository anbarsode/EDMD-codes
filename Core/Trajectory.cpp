// Update documentation
// Test the trajectories wrt other codes
// Trajectories are inaccurate. Must implement RK4 or higher order integration

template <typename T = double>
T rDotSquared(T r, T params[], Potential<T> V)
{
    /*
    Calculates (dr/dt)^2. See Goldstein's book on Classical Mechanics.
    */
    
    T u2 = params[0], ub = params[1], u2b2 = params[2], two_by_m = params[3];
    T KE = u2 - u2b2 / r / r - two_by_m * V.Continuous(r);
    //std::cout << r << " " << KE << std::endl;
    //if(KE < 0) std::cout << "Warning: KE < 0" << std::endl;
    return KE;
}


template <typename T = double>
T BisectionRmin(T params[], T rhi, T dr, Potential<T> V)
{
    /*
    Finds the root of (dr/dt)^2 function using the bisection method
    Resolution is V.dr
    Default lower bracket is 0 initially
    Returns the upper estimate so as to avoid the classically forbidden region
    */
    
    T rlo = 0, rmin;
    
    while(rhi - rlo > dr)
    {
        T fhi = rDotSquared(rhi, params, V);
        T flo = rDotSquared(rlo, params, V);
        if(signbit(fhi) == signbit(flo))
        {
            std::cout << "Please provide a better initial bracket" << std::endl;
            return 0;
        }
        else
        {
            rmin = rhi;
            T rmid = rlo*0.5 + rhi*0.5;
            T fmid = rDotSquared(rmid, params, V);
            if(signbit(fhi) != signbit(fmid)) rlo = rmid;
            else rhi = rmid;
        }
    }
    return rmin;
}


template <typename T = double>
void Simpson13Integrator(T params[], T xmin, T xmax, T dx_min, Potential<T> V, T &Dt, T &Dth)
{
    /*
    Integrates the equations of motion numerically using Simpson's 1/3rd rule
    Radial resolution is at least as fine as dr, and close to it
    */
    
    int Nsteps = floor((xmax - xmin) / dx_min);
    if(Nsteps < 2) Nsteps = 2;
    Nsteps -= Nsteps % 2;
    T dx = (xmax - xmin) / Nsteps;
    
    T dtdr, r;
    
    T dt = 0, dth = 0;
    
    for(int i=0; i<=Nsteps; i++)
    {
        r = xmin + dx * i;
        dtdr = pow(rDotSquared(r, params, V),-0.5);
        if(i == 0 || i == Nsteps)
        {
            dt += dtdr;
            dth += dtdr * params[1] / r / r;
        }
        else if(i % 2 == 0)
        {
            dt += 2 * dtdr;
            dth += 2 * dtdr * params[1] / r / r;
        }
        else
        {
            dt += 4 * dtdr;
            dth += 4 * dtdr * params[1] / r / r;
        }
    }
    
    dt *= dx / 3.0;
    dth *= dx / 3.0;
    
    Dt += fabs(dt);
    Dth += fabs(dth);
} //Inaccurate
    
    
template <typename T = double>
void AdRK4Integrator(T params[], T xmin, T xmax, T dx_min, Potential<T> V, T &Dt, T &Dth)
{
    /*
    Integrates the equations of motion numerically using an adaptive RK4 algorithm
    Radial resolution is at least as fine as dr
    */
    T dtdr, r, dx;
    T c0, c1, c3;
    T d0, d1, d3;
    int maxiter = 8, iter=-1;
    T rtol = 1e-2;
    
    int Nsteps = floor((xmax - xmin) / dx_min);
    if(Nsteps < 2) Nsteps = 2;
    Nsteps -= Nsteps % 2;
    
    T dtc, dthc, dtf=0, dthf=0;
    do
    {
        dtc = dtf;
        dthc = dthf;
        dtf = 0;
        dthf = 0;
        dx = (xmax - xmin) / Nsteps;
        for(int i=0; i<Nsteps; i++)
        {
            r = xmin + dx * i;
            dtdr = pow(rDotSquared(r, params, V),-0.5);
            c0 = dx * dtdr;
            d0 = c0 * params[1] / r / r;
            
            r = xmin + dx * i + dx / 2.0;
            dtdr = pow(rDotSquared(r, params, V),-0.5);
            c1 = dx * dtdr;
            d1 = c1 * params[1] / r / r;
            
            r = xmin + dx * i + dx;
            dtdr = pow(rDotSquared(r, params, V),-0.5);
            c3 = dx * dtdr;
            d3 = c3 * params[1] / r / r;
            
            dtf += (c0 + 4.0 * c1 + c3);
            dthf += (d0 + 4.0 * d1 + d3);
        }
        Nsteps *= 2;
        iter++;
    }while(fabs(dtc/dtf-1) > rtol && fabs(dthc/dthf-1) > rtol && iter < maxiter);
    //if(iter==maxiter) std::cout << "Warning: exceeded maximum number of iterations during adaptive RK4 integration" << std::endl;
    
    dtf /= 6.0;
    dthf /= 6.0;
    
    Dt += fabs(dtf);
    Dth += fabs(dthf);
}




template <class T = double>
class Trajectory
{
    public:
        T u, b, m, dr;
        Potential<T> V; // an object containing cont pot, Nsteps, Rsteps, Vsteps

        T ContRmin, ContDtEnd, ContDthEnd;
        std::vector<T> ContR, ContDt,  ContDth;

        T DiscRmin, DiscDtEnd, DiscDthEnd;
        std::vector<T> DiscR, DiscDt, DiscDth;
        std::vector<T> k1, k2;
        
        std::vector<T> DiffDt, DiffDth;
        std::vector<T> dDtidVi, dDthidVi;
        std::vector<T> dDtidRi, dDtidRim1;
        std::vector<T> dDthidRi, dDthidRim1;


        Trajectory(Potential<T> V_, T u_, T b_, T m_, T dr_)
        {
            V = V_;
            u = u_;
            b = b_;
            m = m_;
            dr = dr_;
            
            ContR.clear();
            ContDt.clear();
            ContDth.clear();
            
            DiscR.clear();
            DiscDt.clear();
            DiscDth.clear();
            k1.clear();
            k2.clear();
            
            DiffDt.clear();
            DiffDth.clear();
            dDtidVi.clear();
            dDthidVi.clear();
            dDtidRi.clear();
            dDtidRim1.clear();
            dDthidRi.clear();
            dDthidRim1.clear();
        }
        
        
        void Clear(void)
        {
            ContR.clear();
            ContDt.clear();
            ContDth.clear();
            
            DiscR.clear();
            DiscDt.clear();
            DiscDth.clear();
            k1.clear();
            k2.clear();
            
            DiffDt.clear();
            DiffDth.clear();
            dDtidVi.clear();
            dDthidVi.clear();
            dDtidRi.clear();
            dDtidRim1.clear();
            dDthidRi.clear();
            dDthidRim1.clear();
        }
        

        void ComputeContTraj(void)
        {
            //Compute (for the first time) cont traj to get ContRmin(maybe a separate function), ContDt, ContDth for ContR
            //Defining separate ContR and DiscR because these R's may be different from Rsteps of V due to different Rmin
            // r_steps must be a decreasing array from Rcut to 0
            // dr needs to be defined
            
            T u2 = u * u, ub = u * b, u2b2 = ub * ub, two_by_m = 2.0 / m;
            T params[4] = {u2, ub, u2b2, two_by_m};
            
            ContRmin = BisectionRmin(params, V.Rsteps[0], dr, V);
            
            //std::cout << "ContRmin = " << ContRmin << std::endl;
            
            int Ntraj = 1;
            for(int i=0; i<V.Nsteps; i++)
                if(ContRmin < V.Rsteps[i]) Ntraj++;
            
            ContR.reserve(Ntraj);
            ContDt.reserve(Ntraj);
            ContDth.reserve(Ntraj);
            
            for(int i=0; i<Ntraj; i++)
            {
                ContR.push_back(0);
                ContDt.push_back(0);
                ContDth.push_back(0);
            }
            
            ContR[0] = V.Rsteps[0]; //Initial r = Rcut
            ContDt[0] = 0; // Initial t = 0 by choice
            ContDth[0] = 0; // Initial theta = 0 by choice
            
            T Dt, Dth;
            
            for(int i=1; i<Ntraj-1; i++)
            {
                Dt = 0;
                Dth = 0;
                AdRK4Integrator<T>(params, V.Rsteps[i-1], V.Rsteps[i], dr, V, Dt, Dth);
                ContR[i] = V.Rsteps[i];
                ContDt[i] = Dt;
                ContDth[i] = Dth;
            }
            
            Dt = 0;
            Dth = 0;
            AdRK4Integrator<T>(params, V.Rsteps[Ntraj-2], ContRmin, dr, V, Dt, Dth);
            ContR[Ntraj-1] = ContRmin;
            ContDt[Ntraj-1] = Dt;
            ContDth[Ntraj-1] = Dth;
        }
        

        void ReComputeContTraj(std::vector<T> NewContR)
        {
            //recalc using already known info from prev calc
            // skipping for now
            ;
        }
        

        void ComputeDiscTraj(void)
        {
            // r_steps must be a decreasing array from Rcut to 0
            // V_steps[i] is the potential between r_steps[i-1] to r_steps[i]
            
            T u2 = u * u, ub = u * b;
            
            DiscR.push_back(V.Rsteps[0]); //Initial r = Rcut
            DiscDt.push_back(0); // Initial t = 0 by choice
            DiscDth.push_back(0); // Initial theta = 0 by choice
            k1.push_back(0);
            k2.push_back(0);
            
            T Dt, Dth, k1i_sq, k1i, k2i, r1_sq, r2_sq;
            
            for(int i=1; i<V.Nsteps; i++)
            {
                Dt = 0;
                Dth = 0;
                
                r1_sq = V.Rsteps[i-1] * V.Rsteps[i-1];
                r2_sq = V.Rsteps[i] * V.Rsteps[i];
                k1i_sq = 1 / (u2 - 2.0 / m * V.Vsteps[i]);
                if(k1i_sq <= 0) break;
                
                k1i = sqrt(k1i_sq);
                k2i = ub * ub * k1i_sq;
                
                k1.push_back(k1i);
                k2.push_back(k2i);
                
                if(r2_sq >= k2i && r1_sq >= k2i)
                {
                    Dt += fabs(k1i * (sqrt(r2_sq - k2i) - sqrt(r1_sq - k2i)));
                    if(b == 0) Dth += 0;
                    else Dth += fabs(atan(sqrt(r2_sq / k2i - 1)) - atan(sqrt(r1_sq / k2i - 1)));
                    DiscR.push_back(V.Rsteps[i]);
                    DiscDt.push_back(Dt);
                    DiscDth.push_back(Dth);
                }
                else if(r1_sq >= k2i)
                {
                    //std::cout << "DiscRmin = " << sqrt(k2i) << std::endl;
                    DiscRmin = sqrt(k2i);
                    Dt += fabs(k1i * sqrt(r1_sq - k2i));
                    if(b == 0) Dth += 0;
                    else Dth += fabs(atan(sqrt(r1_sq / k2i - 1)));
                    DiscR.push_back(sqrt(k2i));
                    DiscDt.push_back(Dt);
                    DiscDth.push_back(Dth);
                    break;
                }
                else break;
            }
        }
        
        
        void ComputeTrajectories(void)
        {
            Clear();
            ComputeContTraj();
            ComputeDiscTraj();
        }
        
        
        void EndPointValues(bool RecomputeTraj = true)
        {
            if(RecomputeTraj == true) ComputeTrajectories();
            
            ContDtEnd = 0;
            ContDthEnd = 0;
            for(int i=0; i<ContR.size(); i++)
            {
                ContDtEnd += ContDt[i];
                ContDthEnd += ContDth[i];
            }
            ContDtEnd *= 2.0;
            ContDthEnd *= 2.0;
            
            DiscDtEnd = 0;
            DiscDthEnd = 0;
            for(int i=0; i<DiscR.size(); i++)
            {
                DiscDtEnd += DiscDt[i];
                DiscDthEnd += DiscDth[i];
            }
            DiscDtEnd *= 2.0;
            DiscDthEnd *= 2.0;
        }
        
        
        void TurningPointCleanup(T BigNum, T SmallNum)
        {
            // It is not necessary that the discrete trajectory will end in the same spherical shell in which the continuous trajectory ends.
            //In that case, there may be trouble when defining gradients and differences later on.
            //This function cleans up the trajectories near their turning points to make them suitable for gradient descent.
            //Note that most of the changes are ad-hoc and may be unphysical in some cases.
            
            if(ContR.size() != 0 && DiscR.size() != 0)
            {
                if(ContRmin != DiscRmin && ContR.end()[-2] != DiscR.end()[-2])
                {
                    if(ContRmin < DiscRmin)
                    {
                        //Cont penetrates deeper
                        //Vi should made smaller, Ri should be made larger (Rs are decr)
                        T k2i, Ri, Rim1;

                        ContR.pop_back();
                        ContDt.pop_back();
                        ContDth.pop_back();

                        DiscR.end()[-1] = ContR[DiscR.size()-1];
                        //No change in Dt, Dth since we have some Disc penetration into this shell, so using those values instead of SmallNum

                        while(DiscR.size() < ContR.size())
                        {
                            DiscR.push-back(ContR[DiscR.size()]);
                            DiscDt.push_back(SmallNum);
                            DiscDth.push_back(SmallNum);
                            k1.push_back(BigNum); //Not sure about this

                            Ri = DiscR.end()[-1];
                            Rim1 = DiscR.end()[-2];
                            k2i = 0.5 * Ri * Ri * (Ri *Ri + 2.0 * Rim1 * Rim1) / (Ri *Ri + Rim1 * Rim1);
                            k2.push_back(k2i);
                        }
                    }
                    else
                    {
                        //Disc penetrates deeper
                        //Vi should be made larger

                        DiscR.pop_back();
                        DiscDt.pop_back();
                        DiscDth.pop_back();

                        ContR.end()[-1] = DiscR[ContR.size()-1];
                        ContDt.end()[-1] = BigNum;
                        ContDth.end()[-1] = BigNum; //Not sure about this

                        while(ContR.size() < DiscR.size())
                        {
                            ContR.push-back(DiscR[ContR.size()]);
                            ContDt.push_back(BigNum);
                            ContDth.push_back(BigNum); //Not sure about this
                        }
                    }
                }
            }
            else std::cout << "Trajectory is empty" << std::endl;
        } //untested as of 3-5-22


        void ComputeDifferences(bool relative = true)
        {
            //find the error between cont and disc trajectories
            
            if(relative == true)
            {
                for(int i=0; i<ContR.size(); i++)
                {
                    DiffDt.push_back(DiscDt[i] / ContDt[i] - 1.0);
                    DiffDth.push_back(DiscDth[i] / ContDth[i] - 1.0);
                }
            }
            else
            {
                for(int i=0; i<ContR.size(); i++)
                {
                    DiffDt.push_back(DiscDt[i] - ContDt[i]);
                    DiffDth.push_back(DiscDth[i] - ContDth[i]);
                }
            }
            
        } //untested as of 3-5-22
        
        
        void ComputePartialDerivatives(void)
        {
            // compute partial derivs of Dt and Dth wrt Vi and Ri
            
            dDtidRi.reserve(DiscR.size());
            dDthidRi.reserve(DiscR.size());
            dDtidRim1.reserve(DiscR.size());
            dDthidRim1.reserve(DiscR.size());
            dDtidVi.reserve(DiscR.size());
            dDthidVi.reserve(DiscR.size());

            dDtidRi.push_back(0);
            dDthidRi.push_back(0);
            dDtidRim1.push_back(0);
            dDthidRim1.push_back(0);
            dDtidVi.push_back(0);
            dDthidVi.push_back(0);

            for(int i=1; i<DiscR.size(); i++)
            {
                T Risq, Rim1sq, k2ByRoot;
                Risq = DiscR[i] * DiscR[i];
                Rim1sq = DiscR[i-1] * DiscR[i-1];
                k2ByRoot = k2[i] / sqrt((Risq - k2[i]) * (Rim1sq - k2[i]));

                dDtidRi.push_back(-DiscR[i] * k1[i] / sqrt(Risq - k2[i]));
                dDtidRim1.push_back(DiscR[i-1] * k1[i] / sqrt(Rim1sq - k2[i]));
                dDthidRi.push_back(u * b * dDtidRi[i] / Risq);
                dDthidRim1.push_back(u * b * dDtidRim1[i] / Rim1sq);
                dDtidVi.push_back(DiscDt[i] * k1[i] * k1[i] / m * (1 - k2ByRoot));
                dDthidVi.push_back(-DiscDt[i] * k2ByRoot / (m * u * b));
            }
        } //untested as of 3-5-22
        
        
        void ComputeTrajDiffDeriv(T BigNum, T SmallNum, bool RelativeDiff = true)
        {
            ComputeTrajectories();
            TurningPointCleanup(BigNum, SmallNum);
            ComputeDifferences(RelativeDiff);
            ComputePartialDerivatives();
        } //untested as of 3-5-22
};
