#include "./PlasProfile.h"

#include <armadillo>

using namespace std;
using namespace arma;

void GenAtomDistribution(vec & Nat, cx_vec const & Z, 
                    double Nat0, double AbsorberWidth, double SourceWidth)
    {
        int Nz = Z.n_elem;
        Nat.zeros(Nz);
        double Zmin = real(Z(0)), Zmax = real(Z(Nz-1)), L = Zmax - Zmin; 
        int minind = floor(double(Nz)*(AbsorberWidth+SourceWidth)/L) , 
        maxind = Nz - floor(double(Nz)*AbsorberWidth/L);
        for (int i = minind; i< maxind; i++)
        {
            //Nat(i) = Nat0;
            if (real(Z(i) - Z(minind) ) < M_PI) 
                Nat(i) = real(  1.0 - cos( Z(i) - Z(minind) )  ) * Nat0 / 2.0;
            if (real( Z(maxind) - Z(i) ) < M_PI) 
                Nat(i) = real(  1.0 - cos( Z(maxind) - Z(i) )  ) * Nat0 / 2.0;
        }

    };