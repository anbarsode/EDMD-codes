template <class T = double>
class Optimizer
{
    public:
        T u, b, m, dr, BigNum=1e-1, SmallNum=1e-6;
        RNG<T> *uDist, *bDist;
        Potential<T> V;
        T Loss;
        
        Optimizer(Potential<T> V_, int uPDF, std::vector<T> uParams, int bPDF, std::vector<T> bParams, T m_ = -1, T dr_ = 1e-6)
        {
            V = V_;
            uDist = new RNG<T>(uPDF, uParams);
            bDist = new RNG<T>(bPDF, bParams);
            m = m_;
            dr = dr_;
        }
        
        T ComputeLoss1(std::mt19937_64& rng, int Neval = 1e3) // EPD without turning point cleanup, relative L2
        {
            T L = 0;
            for(int i=0; i<Neval; i++)
            {
                u = uDist->Dist(rng);
                b = bDist->Dist(rng);
                Trajectory<double> traj(V, u, b, m, dr);
                traj.EndPointValues(true, false);
                if(traj.ContDtEnd != 0) L += pow(traj.DiscDtEnd / traj.ContDtEnd - 1, 2);
                else L += BigNum;
                if(traj.ContDthEnd != 0) L += pow(traj.DiscDthEnd / traj.ContDthEnd - 1, 2);
                else L += BigNum;
            }
            Loss = L / Neval;
            return L / Neval;
        }
        
        T ComputeLoss2(std::mt19937_64& rng, int Neval = 1e3) // EPD with turning point cleanup, relative L2
        {
            T L = 0;
            for(int i=0; i<Neval; i++)
            {
                u = uDist->Dist(rng);
                b = bDist->Dist(rng);
                Trajectory<double> traj(V, u, b, m, dr);
                traj.ComputeTrajectories();
                traj.TurningPointCleanup(BigNum, SmallNum);
                traj.EndPointValues(false, true);
                if(traj.ContDtEnd != 0) L += pow(traj.DiscDtEnd / traj.ContDtEnd - 1, 2);
                else L += BigNum;
                if(traj.ContDthEnd != 0) L += pow(traj.DiscDthEnd / traj.ContDthEnd - 1, 2);
                else L += BigNum;
                //std::cout << L << std::endl;
            }
            Loss = L / Neval;
            return L / Neval;
        }
        
        T ComputeLoss3(std::mt19937_64& rng, int Neval = 1e3) // Traj dev without turning point cleanup, relative L2
        {
            T L = 0;
            for(int i=0; i<Neval; i++)
            {
                u = uDist->Dist(rng);
                b = bDist->Dist(rng);
                Trajectory<double> traj(V, u, b, m, dr);
                traj.ComputeTrajectories();
                traj.DiscardTurningPoint();
                traj.ComputeDifferences(true, BigNum);
                for(int j=0; j<traj.DiffDt.size(); j++) L += pow(traj.DiffDt[j], 2) + pow(traj.DiffDth[j], 2);
            }
            Loss = L / Neval;
            return L / Neval;
        }
        
        T ComputeLoss4(std::mt19937_64& rng, int Neval = 1e3) // Traj dev with turning point cleanup, relative L2
        {
            T L = 0;
            for(int i=0; i<Neval; i++)
            {
                u = uDist->Dist(rng);
                b = bDist->Dist(rng);
                Trajectory<double> traj(V, u, b, m, dr);
                traj.ComputeTrajectories();
                traj.TurningPointCleanup(BigNum, SmallNum);
                traj.ComputeDifferences(true, BigNum);
                for(int j=0; j<traj.DiffDt.size(); j++) L += pow(traj.DiffDt[j], 2) + pow(traj.DiffDth[j], 2);
            }
            Loss = L / Neval;
            return L / Neval;
        }
}; //Completely untested as of 27-5-2022
