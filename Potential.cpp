template <class T>
class Potential
{
    public:
        int Nsteps;
        T Rcut = 3.0; // To be edited at compile time by the user
        T Vcut = 0.0; // TBEACTBTUT
        
        std::vector<T> Rsteps, Vsteps;


        T Continuous(T r) // TBEACTBTU
        {
            /*T k = 1.0; // TBEACTBTU
            return k/r; // maybe only when r <= Rcut, Vcut otherwise?*/
            
            if(r < Rcut)
            {
                T V;
                T A = 1.0, B = 1.0; // TBEACTBTUT
                
                T r6 = pow(r, -6);
                V = (A * r6 - B) * r6; // Lennard Jones
                
                return V - Vcut;
            }
            else return 0.0;
        }


        void Read_Discrete(std::string fname)
        {
            //Read Nsteps, Vsteps, Rsteps from ascii file
            ;
        }


        void Write_Discrete(std::string fname)
        {
            //Write Nsteps, Vsteps, Rsteps from ascii file
            ;
        }


        T Discrete(T r)
        {
            //Return disc pot value based on Rsteps and Vsteps
            return 0;
        }
        
        
        
};
