#ifndef PLASPROFILE_H // include guard
#define PLASPROFILE_H

#include <armadillo>
#include <complex>

using namespace std;
using namespace arma;

struct MyPerturbation
{
    int Hx;
    int Hz;
    int Nz;
    int Nx;
    double Amplitude;
    Mat<double> Shape;
    Col<double> ZShape;
    MyPerturbation(int iNx, int iNz, int iHx, int iHz, double iAmplitude, cx_rowvec x, cx_colvec z, double Lf)
    {
        Hx = iHx;
        Hz = iHz;
        Nz = iNz;
        Nx = iNx;
        Amplitude = iAmplitude;
        Shape.zeros(Nz, Nx);

        Mat<double> pert_z_Shape,
            pert_x_Shape;

        pert_x_Shape.zeros(Nz, Nx);
        pert_z_Shape.zeros(Nz, Nx);

        // Shape of the z-enevelope;
        ZShape.zeros(Nz);
        for (int i = 0; i < Nz; i++)
        {
            if (abs(z(i)) <= Lf)
                ZShape(i) = real(sin(M_PI * (z(i) - Lf) / 2.0 / Lf) * sin(M_PI * (z(i) - Lf) / 2.0 / Lf));
        }

        // X-shape of perturbation;
        cx_double dx = x(1) - x(0), Lx = (x(Nx - 1) - x(0) + dx);
        if (Hx == 0 && abs(Hz) > 0)
        {
            pert_x_Shape.ones(Nz, Nx);
        }
        if (Hx <= -1)
        {
            Mat<double> Rand;
            arma_rng::set_seed_random();
            Rand.randu(1, 1);
            pert_x_Shape.each_row() = -real(cos(2.0 * M_PI * Rand(0) + Hx * 2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx));
        }
        if (Hx == 1)
        {
            pert_x_Shape.each_row() = -real(cos(2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx));
        }
        if (Hx > 1 && Hx < 100)
        {
            Mat<double> Rand;
            arma_rng::set_seed_random();
            Rand.randu(1, Hx);
            for (int j = 0; j < Hx; j++)
            {
                pert_x_Shape.each_row() += -real(cos(2.0 * j * M_PI * (x + dx / 2.0 - x(0)) / Lx + 2.0 * M_PI * Rand(j)));
                //cout << Rand(j) << endl;
            }
        }
        if (Hx == 100)
        {
            pert_x_Shape.each_row() = pow(real(0.5 - 0.5 * cos(2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx)), 100);
        }
        if (Hx == 101)
        {
            Mat<double> Rand;
            arma_rng::set_seed_random();
            Rand.randu(1, 1);
            pert_x_Shape.each_row() = pow(real(0.5 - 0.5 * cos(2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx + 2.0 * M_PI * Rand(0))), 100);
        }
        // Z-shape of perturbation;
        cx_double dz = z(1) - z(0);

        if (Hz == 0)
        {
            pert_z_Shape.ones(Nz, Nx);
        }
        if (Hz <= -1)
        {
            Mat<double> Rand;
            arma_rng::set_seed_random();
            Rand.randu(1, 1);
            pert_z_Shape.each_col() = -real(cos(2.0 * M_PI * Rand(0) + Hz * 2.0 * M_PI * (z + dz / 2.0) / Lx));
        }
        if (Hz == 1)
        {
            pert_z_Shape.each_col() = -real(cos(2.0 * M_PI * (z + dz / 2.0) / Lx));
        }
        if (Hz > 1 && Hz < 100)
        {
            Mat<double> Rand;
            arma_rng::set_seed_random();
            Rand.randu(1, Hz);
            for (int j = 0; j < Hz; j++)
            {
                pert_z_Shape.each_col() += -real(cos(2.0 * j * M_PI * (z + dz / 2.0) / Lx + 2.0 * M_PI * Rand(j)));
            }
        }
        if (Hz == 100)
        {
            int Nz1 = 0, Nz2 = Nz - 1;

            while (abs(z(Nz1)) > abs(Lx/2.0))
                Nz1++;

            while (abs(z(Nz2)) > abs(Lx/2.0))
                Nz2--;

            for (Nz1; Nz1 <= Nz2; Nz1++) {             
                double val = pow(real(0.5 + 0.5 * cos(2.0 * M_PI * (z(Nz1)) / Lx)), 100);
                pert_z_Shape.row(Nz1).for_each([val](double &V) { V = val; }); } //
        }

        Shape = pert_z_Shape % pert_x_Shape;
        if (Shape.max() > 0.0)
            Shape = Shape / Shape.max() * Amplitude;
        for (int i = 0; i < Nx; i++)
            Shape.col(i) = Shape.col(i) % ZShape;
    }
};

void GenAtomDistribution(cx_vec &Nat, cx_vec const &Z,
                         double Nat0, double AbsorberWidth, double SourceWidth);

void GenTstConcentration(mat &Co, cx_colvec const &Z, cx_rowvec const &X);
void GenTstConcentration2(mat &Co, cx_colvec const &Z, cx_rowvec const &X);
void GenTstConcentrationSharp(mat &Co, cx_colvec const &Z, cx_rowvec const &X);
void GenTstConcentrationHomX(mat &Co, cx_colvec const &Z, cx_rowvec const &X);
void GenTstConcentrationHomXSharp(mat &Co, cx_colvec const &Z, cx_rowvec const &X);

#endif