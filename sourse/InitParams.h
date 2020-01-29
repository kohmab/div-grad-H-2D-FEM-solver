#ifndef INITPARAMS_H // include guard
#define INITPARAMS_H

#include "./MyPhys.h"
#include <cmath>

struct _Silica
{
    double Permittivity = 2.1;
    double AtomConcentration = 2.1e22;  /* [cm-3] */
    double RecombinationTime = 150e-15; /* [s] */
    double SigmaInside = 1.5e-28;            /* [cm^3/W^3/s] */
    double SigmaOutside = 1.0e-69;            /* [cm^9/W^6/s] */
    double BandGapInside = 5.2 * 1.6e-12;     /* [erg] */
    double BandGapOutside = 9.1 * 1.6e-12;     /* [erg] */
    int PhotonsNumInside = 3;
    int PhotonsNumOutside = 6;
    double DiffusuinCoofitient = 10.0;  /* [cm2/s] */
} Silica;

struct _Field
{
    double Frequency = 2.35e15;                                                                                          /* [s-1] */
    double MaximumIntensity = 1e13;                                                                                      /* [W/cm2] */
    double CriticalConcentration = Electron.Mass * Frequency * Frequency / 4.0 / pi / Electron.Charge / Electron.Charge; /* [cm-3] */
    double WaveLength = 2 * pi * LightSpeed / Frequency / sqrt(Silica.Permittivity);                                     /* [cm] */
    double WaveNumber = 2 * pi / WaveLength;
} Field;

#endif /* INITPARAMS_H */