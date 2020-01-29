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
void GenTstConcentration(mat & Co, cx_colvec const &  Z, cx_rowvec const & X)
{
    cx_rowvec R;
    int Nx = X.n_elem, Nz = Z.n_elem;
    Co.zeros(Nz,Nx);
    R.ones(Nx);
    double CoLz = 1.5, CoMa = 0.333, CoLx = real(X(Nx-1))/6.0;
    for (int i = 0; i < Nz; i++)
    {
        if (abs(Z(i)) < CoLz && abs(Z(i)) > CoLz - 2.0*M_PI)
        {
            double arg = ( abs(Z(i)) - CoLz + 2.0*M_PI )/2.0;
            Co.row(i) = real(R * (1.0 + cos(arg))/2.0 * CoMa);
        }

        if (abs(Z(i)) < CoLz - 2.0*M_PI)
        {
            Co.row(i) = real(R * CoMa);
        }
    }

    for (int i = 0; i < Nx; i++)
    {
        if (abs(X(i)) < CoLx && abs(X(i)) > CoLx/2.0 )
        {
            double arg = (CoLx - abs(X(i)))/(CoLx)*M_PI;
            Co.col(i) = real(Co.col(i) * sin(arg)*sin(arg));
        }

        if (abs(X(i)) < CoLx/2.0)
        {
            Co.col(i) = Co.col(i);
        }

        if (abs(X(i)) >= CoLx)
        {
           Co.col(i) *= 0.0;
        }
    }   
}

void GenTstConcentration2(mat & Co, cx_colvec const &  Z, cx_rowvec const & X)
{
    cx_rowvec R;
    int Nx = X.n_elem, Nz = Z.n_elem;
    Co.zeros(Nz,Nx);
    R.ones(Nx);
    double CoLz = 1.5, CoMa = 0.333, CoLx = real(X(Nx-1))/5.0;
    for (int i = 0; i < Nz; i++)
    {
        if (abs(Z(i)) < CoLz )
        {
            double arg = abs(Z(i))*M_PI / CoLz;
            Co.row(i) = real(R * (1.0 + cos(arg))/2.0 * CoMa);
        }        
    }

    for (int i = 0; i < Nx; i++)
    {
        if (abs(X(i) - (X(1)-X(0))/2.0) < CoLx )
        {
            double arg = abs(X(i) - (X(1)-X(0))/2.0 ) / CoLx *M_PI;
            Co.col(i) = real(Co.col(i)) * (1.0 + cos(arg))/2.0;
        }

       

        if (abs(X(i) - (X(1)-X(0))/2.0) >= CoLx)
        {
           Co.col(i) *= 0.0;
        }
    }   
}

void GenTstConcentrationSharp(mat & Co, cx_colvec const &  Z, cx_rowvec const & X)
{
    cx_rowvec R;
    int Nx = X.n_elem, Nz = Z.n_elem;
    Co.zeros(Nz,Nx);
    R.ones(Nx);
    double CoLz = 2.5, CoMa = 4.0, CoLx = real(X(Nx-1))/2.0;
    for (int i = 0; i < Nz; i++)
    {
        if (abs(Z(i)) < CoLz)
        {
            Co.row(i) = real(R * CoMa);
        }
    }

    for (int i = 0; i < Nx; i++)
    {
        if (abs(X(i)) < CoLx/3.0)
        {
            Co.col(i) = Co.col(i);
        }

        if (abs(X(i)) >= CoLx/3.0)
        {
           Co.col(i) *= 0.0;
        }
    }   
}

void GenTstConcentrationHomX(mat & Co, cx_colvec const &  Z, cx_rowvec const & X)
{
    cx_rowvec R;
    int Nx = X.n_elem, Nz = Z.n_elem;
    Co.zeros(Nz,Nx);
    R.ones(Nx);
    double CoLz = 2.5*M_PI, CoMa = 0.1;
    for (int i = 0; i < Nz; i++)
    {
        if (abs(Z(i)) < CoLz && abs(Z(i)) > CoLz - 2.0*M_PI)
        {
            double arg = ( abs(Z(i)) - CoLz + 2.0*M_PI )/2.0;
            Co.row(i) = real(R * (1.0 + cos(arg))/2.0 * CoMa);
        }

        if (abs(Z(i)) < CoLz - 2.0*M_PI)
        {
            Co.row(i) = real(R * CoMa);
        }
    }    
}

void GenTstConcentrationHomXSharp(mat & Co, cx_colvec const &  Z, cx_rowvec const & X)
{
    cx_rowvec R;
    int Nx = X.n_elem, Nz = Z.n_elem;
    Co.zeros(Nz,Nx);
    R.ones(Nx);
    double CoLz = 2.5*M_PI, CoMa = 0.3;
    for (int i = 0; i < Nz; i++)
    {
        if (abs(Z(i)) < CoLz - 2.0*M_PI)
        {
            Co.row(i) = real(R * CoMa);
        }
    }    
}

    //Plot_ABS2D(gpN, Co, "N"); */