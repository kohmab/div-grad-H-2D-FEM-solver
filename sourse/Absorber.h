#pragma once
#ifndef ABSORBER_H // include guard
#define ABSORBER_H
#include <armadillo>
#include <complex>
#include "./MyPhys.h"

using namespace std;
using namespace arma;

struct MyAbsorber
{
    complex<double> A, B;
    double Width;
    cx_vec ShapeN;
    cx_vec ShapeH;
    //cx_vec Shape;
    MyAbsorber(cx_vec Zn, cx_vec Zh, double W, complex<double> Ai, complex<double> Bi)
    {
        //      complex<double> I(0,1);
        A = Ai;
        B = Bi;
        Width = W;
        int Nz = Zn.n_elem;
        ShapeN.zeros(Nz);
        for (int i = 0; i < Nz; i++)
        {
            if (real(Zn(i) - Zn(0)) < real(W))
            {
                complex<double> X = (Zn(i) - Zn(0)) * pi / W;
                ShapeN(i) = (1.0 - cos(X) * cos(X)) * A;
            }
            if (real(Zn(i) - Zn(0)) < real(W / 2.0))
            {
                ShapeN(i) = A;
            }

            if (real(Zn(Nz - 1) - Zn(i)) < real(W))
            {
                complex<double> X = (Zn(Nz - 1) - Zn(i)) * pi / W;
                ShapeN(i) = (1.0 - cos(X) * cos(X)) * A;
            }
            if (real(Zn(Nz - 1) - Zn(i)) < real(W / 2.0))
            {
                ShapeN(i) = A;
            }
        }
        Nz = Zh.n_elem;
        ShapeH.zeros(Nz);
        for (int i = 0; i < Nz; i++)
        {
            if (real(Zh(i) - Zh(0)) < real(W))
            {
                complex<double> X = (Zh(i) - Zh(0)) * pi / W;
                ShapeH(i) = (1.0 - cos(X) * cos(X)) * B;
            }
            if (real(Zh(i) - Zh(0)) < real(W / 2.0))
            {
                ShapeH(i) = B;
            }
            if (real(Zn(Nz - 1) - Zn(i)) < real(W))
            {
                complex<double> X = (Zh(Nz - 1) - Zh(i)) * pi / W;
                ShapeH(i) = (1.0 - cos(X) * cos(X)) * B;
            }
            if (real(Zn(Nz - 1) - Zn(i)) < real(W / 2.0))
            {
                ShapeH(i) = B;
            }
        }
    };
};

#endif //ABSORBER_H