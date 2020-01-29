#pragma once
#ifndef MYPHYS_H // include guard
#define MYPHYS_H
#include <complex>

struct _Electron
{
    double Mass = 0.910938356e-27;  /* [Gramm] */
    double Charge = 4.80320425e-10; /* [CGSE] */
} Electron;

double LightSpeed = 2.99792458e10; /* [cm/s] */

double pi = 3.14159265359;

std::complex<double> optp(0.5 / pi, 0.0), I(0.0, 1.0);

#endif /* MYPHYS_H */