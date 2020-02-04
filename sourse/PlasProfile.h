#ifndef PLASPROFILE_H // include guard
#define PLASPROFILE_H

#include <armadillo>
#include <complex>

using namespace std;
using namespace arma;

struct MyPerturbation
{
    int Nh;
    int Nw;
    int Nz;
    int Nx;
    double Amplitude;
    Mat<double> Shape;
    Col<double> ZShape;
    MyPerturbation(int iNx, int iNz, int iNh, int iNw, double iAmplitude, cx_rowvec x, cx_colvec z, double Lf)
    {
        Nw = iNw;
        Nh = iNh;
        Nz = iNz;
        Nx = iNx;
        Amplitude = iAmplitude;
        Shape.zeros(Nz, Nx);

        // Shape of the z-enevelope;
        ZShape.zeros(Nz);
        for (int i = 0; i < Nz; i++)
        {
            if (abs(z(i)) <= Lf)
                ZShape(i) = real(sin(M_PI * (z(i) - Lf) / 2.0 / Lf) * sin(M_PI * (z(i) - Lf) / 2.0 / Lf));
        }

        cx_double dx = x(1) - x(0), Lx = (x(Nx - 1) - x(0) + dx);

        if (Nw == 0)
        {
            if (Nh == 0)
            {
                Shape.ones(Nz, Nx);
            }
            if (Nh <= -1)
            {
                Mat<double> Rand;
                arma_rng::set_seed_random();
                Rand.randu(1, 1);
                Shape.each_row() = -real(cos(2.0 * M_PI * Rand(0) + Nh * 2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx));
            }

            if (Nh == 1)
            {
                Shape.each_row() = -real(cos(2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx));
            }
            if (Nh > 1)
            {
                Mat<double> Rand;
                arma_rng::set_seed_random();
                Rand.randu(1, Nh);
                for (int j = 0; j < Nh; j++)
                {
                    Shape.each_row() += -real(cos(2.0 * j * M_PI * (x + dx / 2.0 - x(0)) / Lx + 2.0 * M_PI * Rand(j)));
                    //cout << Rand(j) << endl;
                }
            }
            if (Nh == 100)
            {
                Shape.each_row() = pow(real(0.5 - 0.5 * cos(2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx)), 100);
            }
            if (Nh == 101)
            {
                Mat<double> Rand;
                arma_rng::set_seed_random();
                Rand.randu(1, 1);
                Shape.each_row() = pow(real(0.5 - 0.5 * cos(2.0 * M_PI * (x + dx / 2.0 - x(0)) / Lx + 2.0 * M_PI * Rand(0))), 100);
            }
        }
        else
        {
            mat X2D(Nz, Nx, fill::zeros), Z2D(Nz, Nx, fill::zeros);
            X2D.each_row() = real(x);
            Z2D.each_col() = real(z);
            for (int h = 1; h <= Nh; h++)
            {
                vec phase(Nw, fill::randu), angle(Nw, fill::randu);
                phase *= 2.0 * M_PI;
                angle *= M_PI;
                double k = 2 * M_PI / real(Lx) * (double)h;
                for (int w = 0; w < Nw; w++)
                {
                    Shape += cos(k * (X2D * cos(angle(w)) + Z2D * sin(angle(w))) + phase(w));
                }
            }
        }

        if (Shape.max() > 0.0)
            Shape = Shape / Shape.max() * Amplitude;
        for (int i = 0; i < Nx; i++)
            Shape.col(i) = Shape.col(i) % ZShape;
    }
};

void GenAtomDistribution(cx_vec &Nat, cx_vec const &Z,
                         double Nat0, double AbsorberWidth, double SourceWidth);

#endif