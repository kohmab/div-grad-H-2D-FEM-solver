#pragma once
#ifndef EFFECTIVEMEDIA_H // include guard
#define EFFECTIVEMEDIA_H
#include <armadillo>
#include <omp.h>
#include <complex>

using namespace std;
using namespace arma;

struct MyMedia
{
    cx_mat eps_x, eps_z, al_x, al_z;
    mat apb, f;

    MyMedia(int Numel,
            double apb0, double f0, double Nc,
            cx_double De,
            double D, double SigmaAOut, mat &N0)
        : SigmaA{SigmaAOut}, N{Numel}, De{De}, Nc{Nc}, Diffn{D}
    {
        a.ones(N, 1);
        b.ones(N, 1);
        b /= apb0;
        eps_x.ones(N, 1);
        eps_z.ones(N, 1);
        al_x.ones(N, 1);
        al_z.ones(N, 1);
        f.ones(N, 1);
        f *= f0;
        apb.ones(N, 1);
        apb *= apb0;
        cx_double eks0 = eks_calc(apb0);
        eks.ones(N, 1);
        eks *= eks0;
        dep_x.ones(N, 1);
        dep_x *= dep_x_calc(eks0);
        dep_z = (1.0 - dep_x) / 2.0;
        for (uword i = 0; i < N; i++)
        {
            eps_x[i] = eps_x_calc(N0[i] * (1.0 + De), i);
            eps_z[i] = eps_z_calc(N0[i] * (1.0 + De), i);
            al_x[i] = al_x_calc(N0[i] * (1.0 + De), i);
            al_z[i] = al_z_calc(N0[i] * (1.0 + De), i);
        }
    }

    void update(mat const &N0, mat const &N1, cx_mat const &Exa, cx_mat const &Eza, double dT)
    {
#pragma omp paralel for
        for (int i = 0; i < N; i++)
        {
            if ((N1[i] > N0[i] + dT * 1.0e-5) && (f[i] < 1.0))
            {
                cx_double eps = 1.0 - (N0[i] + N1[i]) / 2.0 * (1.0 + De);
                double NuInside = 2.0 * (N1[i] - N0[i]) / (N1[i] + N0[i]),
                       NuPole = SigmaA * (pow(abs(Exa[i] * al_x[i] * eps), 2) +
                                          pow(abs(Eza[i] * al_z[i]), 2)),
                       NuEquator = SigmaA * (pow(abs(Exa[i] * al_x[i]), 2) +
                                             pow(abs(Eza[i] * al_z[i] * eps), 2));
                a[i] += sqrt((NuPole + NuInside) / 2.0 * Diffn) * dT;
                b[i] += sqrt((NuEquator + NuInside) / 2.0 * Diffn) * dT;
                f[i] = a[i] * b[i] * b[i] * Nc;
                apb[i] = a[i] / b[i];
                eks[i] = eks_calc(apb[i]);
                dep_x[i] = dep_x_calc(eks[i]);
                dep_z[i] = (1.0 - dep_x[i]) / 2.0;
            }
            if (f[i] > 1.0)
                f[i] = 1.0;
            al_x[i] = al_x_calc(N0[i] * (1.0 + De), i);
            al_z[i] = al_z_calc(N0[i] * (1.0 + De), i);
            eps_x[i] = eps_x_calc(N1[i] * (1.0 + De), i);
            eps_z[i] = eps_z_calc(N1[i] * (1.0 + De), i);
        }
    }

private:
    Mat<double> dep_x, dep_z, a, b;
    cx_mat eks;
    double SigmaA, Nc, Diffn;
    cx_double De;
    int N;
    cx_double inline eps_z_calc(cx_double N, int n)
    {
        return (1.0 - (f[n] + dep_z[n] * (1.0 - f[n])) * N) / (1.0 - dep_z[n] * (1.0 - f[n]) * N);
    }
    cx_double inline eps_x_calc(cx_double N, int n)
    {
        return (1.0 - (f[n] + dep_x[n] * (1.0 - f[n])) * N) / (1.0 - dep_x[n] * (1.0 - f[n]) * N);
    }
    cx_double inline al_z_calc(cx_double N, int n)
    {
        return 1.0 / (1.0 - dep_z[n] * (1.0 - f[n]) * N);
    }
    cx_double inline al_x_calc(cx_double N, int n)
    {
        return 1.0 / (1.0 - dep_x[n] * (1.0 - f[n]) * N);
    }
    double inline dep_x_calc(cx_double eks)
    {
        return (eks == 0.0)
                   ? 1.0 / 3.0
                   : real((1.0 - eks * eks) / eks / eks / eks * (atanh(eks) - eks));
    }
    cx_double inline eks_calc(double apb)
    {
        return sqrt(complex<double>(1.0 - 1.0 / apb / apb));
    }
};

#endif //EFFECTIVEMEDIA_H