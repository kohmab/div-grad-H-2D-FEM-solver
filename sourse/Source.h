#pragma once
#ifndef SOURCE_H // include guard
#define SOURCE_H

#include <armadillo>
#include "./Absorber.h"

using namespace std;
using namespace arma;

struct MySource
{
public:
    cx_mat Shape;
    double Amplitude;
    double ZWidth;
    double ZShift;
    double TWidth;
    double TShift;
    double TimeFunction(double T)
    {
        double d = (tmpf(T) - tmpf0) * tmpfnorm;
        return d < 0.0 ? 0.0 : d;
    };
    double TimeFunctionIntegral(double T, double dT)
    {
        double d = TimeFunction(T) * (dT + (TShift - T) / 2.0 / TWidth / TWidth * dT * dT);
        return d < 0.0 ? 0.0 : d;
    };

    MySource(int type, cx_vec Z, cx_rowvec X, double LF, double ZS, double TW, double TS, double FieldCalibrationConst)
    {
        complex<double> I(0.0, 1.0);
        double Zc;
        double Kx = real(2.0 * M_PI / (X(X.n_elem - 1) - 2.0*X(0)+X(1))) * double(type - 1),
               Kz = real(sqrt(4.0 * M_PI * M_PI - Kx * Kx));
        double ZW = 2.0 * M_PI / Kz;
        Zc = (ZS + ZW + real(Z(0, 0))) / LF;
        Amplitude = 1.0;
        ZWidth = ZW;
        ZShift = ZS;
        TWidth = TW;
        TShift = TS;
        tmpf0 = tmpf(0);
        tmpfnorm = 1.0 / (1.0 - tmpf0);
        Shape.zeros(Z.n_elem, X.n_elem);

        for (int i = 0; i < Z.n_elem; i++)
            for (int j = 0; j < X.n_elem; j++)

            {
                Zc = real(Z(i) - Z(0));
                if ((Zc > ZS) && (Zc < ZS + ZW))
                {
                    /* Shape(i,j) = sin((Z(i) - ZS) / ZW * M_PI) * sin((Z(i) - ZS) / ZW * M_PI) * 2.0 / M_PI *
                               exp(-2.0 * M_PI * I * Z(i)) /
                               FieldCalibrationConst; //0.306498100;  */
                    Shape(i, j) = exp(-Kz * I * Z(i)) / M_PI * cos(Kx * X(j)) /
                                  FieldCalibrationConst; //0.306498100;
                }
            }
    };

private:
    double tmpf(double T)
    {
        return exp(-(T - TShift) * (T - TShift) / TWidth / TWidth / 2.0);
    };
    double tmpf0, tmpfnorm;
};

void CalcSourceAmplitude(MySource &Source, MyAbsorber const &Absorber, cx_vec &TubeFunc, double dZ)
{
    int Nz = TubeFunc.n_elem;
    cx_vec E = zeros<cx_vec>(Nz);
    cx_vec C = ones<cx_vec>(Nz);
    cx_vec F = zeros<cx_vec>(Nz);
    cx_vec al = zeros<cx_vec>(Nz);
    cx_vec bt = zeros<cx_vec>(Nz);

    F.rows(1, Nz - 2) = dZ * dZ * 4.0 * pi * pi * (Source.Shape.rows(2, Nz - 1));
    C.rows(1, Nz - 2) = 2.0 - (1.0 - Absorber.ShapeH.rows(1, Nz - 2)) * dZ * dZ * 4.0 * pi * pi;

    for (int i = 0; i < Nz - 2; i++)
    {
        al(i + 1) = (1.0 + TubeFunc(i) / 2.0) / (C(i) - (1.0 - TubeFunc(i) / 2.0) * al(i));
        bt(i + 1) = (F(i) + (1.0 - TubeFunc(i) / 2.0) * bt(i)) / (C(i) - (1.0 - TubeFunc(i) / 2.0) * al(i));
    }
    E(Nz - 1) = 0.0;
    for (int i = Nz - 1; i > 1; i--)
    {
        E(i - 1) = al(i) * E(i) + bt(i);
    }
    E(0) = 0.0;
    cx_vec Z = linspace<cx_vec>(0, Nz - 1, Nz) * dZ;

    Source.Amplitude = Source.Amplitude;            //max(abs(E));
    Source.Shape = Source.Shape * Source.Amplitude; //0.93696445398592;
}
#endif