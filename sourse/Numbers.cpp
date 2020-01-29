#include <iostream>
#include <armadillo>
//#include <blas.h>
#include <omp.h>

#include "./MyPhys.h"
#include "./SetConfig.h"

using namespace arma;
using arma::Mat;
using arma::span;
using arma::zeros;
using std::cin;
using std::complex;
using std::cout;

int main()
{
    ifstream ConfigFile;
    string ConfigString;

    try
    {
        ConfigFile.open("./config.m");
        throw ConfigFile.is_open();
    }
    catch (bool Op)
    {
        if (Op == false)
            ConfigFile.open("../config.m");
    }

    int nthreads = 4, Nz, Nx;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, nthreads);
    omp_set_num_threads(nthreads);

    for (int i = 0; i < 2; i++)
        getline(ConfigFile, ConfigString);
    GetParam(ConfigString, Nz);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, Nx);

    double ZminLam, ZmaxLam;
    double XminLam, XmaxLam;
    for (int i = 0; i < 6; i++)
        getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ZminLam);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ZmaxLam);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, XminLam);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, XmaxLam);
    double Zmin = ZminLam, Zmax = ZmaxLam, ZLength = Zmax - Zmin, dZ = ZLength / double(Nz + 1);
    double XLength = (XmaxLam - XminLam), dX = XLength / double(Nx),
           Xmin = XminLam, Xmax = XmaxLam - dX;

    double L;
    for (int i = 0; i < 6; i++)
        getline(ConfigFile, ConfigString);
    GetParam(ConfigString, L);
    L = L / 2.0 / M_PI;
    ConfigFile.close();

    //-----------------------------------------------------------------------------------------------------------------------
    int Nmax = Nx * Nz;
    Mat<int> NumbersA, NumbersEps, NeibourA, NeibourEps;
    NumbersA.zeros(Nmax, 2);
    NumbersEps.zeros(Nmax + Nx, 2);
    NeibourA.zeros(Nmax, 9);
    NeibourEps.zeros(Nmax, 16);

    int ind;
    for (int xi = 0; xi < Nx; xi++)
        for (int zi = 0; zi < Nz; zi++)
        {
            ind = zi * Nx + xi;
            NumbersA(ind, 0) = zi;
            NumbersA(ind, 1) = xi;

            NumbersEps(ind, 0) = zi;
            NumbersEps(ind, 1) = xi;

            NeibourA(ind, 0) = (zi == 0) ? -1 : (xi % Nx == 0) ? ind - 1 : ind - Nx - 1;
            NeibourA(ind, 1) = (zi == 0) ? -1 : ind - Nx;
            NeibourA(ind, 2) = (zi == 0) ? -1 : (xi % Nx == Nx - 1) ? ind - 2 * Nx + 1 : ind - Nx + 1;
            NeibourA(ind, 3) = (xi % Nx == 0) ? ind + Nx - 1 : ind - 1;
            NeibourA(ind, 4) = ind;
            NeibourA(ind, 5) = (xi % Nx == Nx - 1) ? ind - Nx + 1 : ind + 1;
            NeibourA(ind, 6) = (zi == Nz - 1) ? -1 : (xi % Nx == 0) ? ind + 2 * Nx - 1 : ind + Nx - 1;
            NeibourA(ind, 7) = (zi == Nz - 1) ? -1 : ind + Nx;
            NeibourA(ind, 8) = (zi == Nz - 1) ? -1 : (xi % Nx == Nx - 1) ? ind + 1 : ind + Nx + 1;

            NeibourEps(ind, 0) = (zi == 0) ? -1 : (xi % Nx <= 1) ? ind - 2 : ind - Nx - 2;
            NeibourEps(ind, 1) = (zi == 0) ? -1 : (xi % Nx == 0) ? ind - 1 : ind - Nx - 1;
            NeibourEps(ind, 2) = (zi == 0) ? -1 : ind - Nx;
            NeibourEps(ind, 3) = (zi == 0) ? -1 : (xi % Nx == Nx - 1) ? ind - 2 * Nx + 1 : ind - Nx + 1;
            NeibourEps(ind, 4) = (xi % Nx <= 1) ? ind + Nx - 2 : ind - 2;
            NeibourEps(ind, 5) = (xi % Nx == 0) ? ind + Nx - 1 : ind - 1;
            NeibourEps(ind, 6) = ind;
            NeibourEps(ind, 7) = (xi % Nx == Nx - 1) ? ind - Nx + 1 : ind + 1;
            NeibourEps(ind, 8) = (xi % Nx <= 1) ? ind + 2 * Nx - 2 : ind + Nx - 2;
            NeibourEps(ind, 9) = (xi % Nx == 0) ? ind + 2 * Nx - 1 : ind + Nx - 1;
            NeibourEps(ind, 10) = ind + Nx;
            NeibourEps(ind, 11) = (xi % Nx == Nx - 1) ? ind + 1 : ind + Nx + 1;
            NeibourEps(ind, 12) = (zi == Nz - 1) ? -1 : (xi % Nx <= 1) ? ind + 3 * Nx - 2 : ind + 2 * Nx - 2;
            NeibourEps(ind, 13) = (zi == Nz - 1) ? -1 : (xi % Nx == 0) ? ind + 3 * Nx - 1 : ind + 2 * Nx - 1;
            NeibourEps(ind, 14) = (zi == Nz - 1) ? -1 : ind + 2 * Nx;
            NeibourEps(ind, 15) = (zi == Nz - 1) ? -1 : (xi % Nx == Nx - 1) ? ind + Nx + 1 : ind + 2 * Nx + 1;
        }
    for (int xi = 0; xi < Nx; xi++)
    {
        ind = Nz * Nx + xi;
        NumbersEps(ind, 0) = Nz;
        NumbersEps(ind, 1) = xi;
    }

    if (Nmax < 100)
    {
        NumbersA.print("Elemetns indexes: ");
        NumbersEps.print("Permittivity indexes: ");

        NeibourA.print("Element neibours: ");
        NeibourEps.print("Permittivity around element: ");
    }
    NumbersA.save(hdf5_name("../Matrices.h5", "/ElementsIndexes"));
    NumbersEps.save(hdf5_name("../Matrices.h5", "/PermittivityIndexes", hdf5_opts::append));
    NeibourA.save(hdf5_name("../Matrices.h5", "/NeighboringElements", hdf5_opts::append));
    NeibourEps.save(hdf5_name("../Matrices.h5", "/PermittivityAround", hdf5_opts::append));

    //-----------------------------------------------------------------------------------------------------------------------
    Mat<double> EpsilonMaskX, EpsilonMaskZ;
    double xpz = dX / dZ, zpx = dZ / dX;
    EpsilonMaskZ.zeros(9, 16);
    EpsilonMaskX.zeros(9, 16);

    Mat<double> Tmp_Row;

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = -3.0 * zpx;
    Tmp_Row(2) = Tmp_Row(0);
    Tmp_Row(8) = Tmp_Row(0);
    Tmp_Row(10) = Tmp_Row(0);
    Tmp_Row(1) = -18.0 * zpx;
    Tmp_Row(9) = Tmp_Row(1);
    Tmp_Row(4) = -26.0 * zpx;
    Tmp_Row(6) = Tmp_Row(4);
    Tmp_Row(5) = -156.0 * zpx;
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskZ.row(0).cols(0, 10) = Tmp_Row.cols(0, 10);
    EpsilonMaskZ.row(2).cols(1, 11) = Tmp_Row.cols(0, 10);
    EpsilonMaskZ.row(6).cols(4, 14) = Tmp_Row.cols(0, 10);
    EpsilonMaskZ.row(8).cols(5, 15) = Tmp_Row.cols(0, 10);

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = 3.0 * zpx;
    Tmp_Row(3) = Tmp_Row(0);
    Tmp_Row(8) = Tmp_Row(0);
    Tmp_Row(11) = Tmp_Row(0);
    Tmp_Row(1) = 21.0 * zpx;
    Tmp_Row(2) = Tmp_Row(1);
    Tmp_Row(9) = Tmp_Row(1);
    Tmp_Row(10) = Tmp_Row(1);
    Tmp_Row(4) = 26.0 * zpx;
    Tmp_Row(7) = Tmp_Row(4);
    Tmp_Row(5) = 182.0 * zpx;
    Tmp_Row(6) = Tmp_Row(5);
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskZ.row(1).cols(0, 11) = Tmp_Row.cols(0, 11);
    EpsilonMaskZ.row(7).cols(4, 15) = Tmp_Row.cols(0, 11);

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = -1.0 * zpx;
    Tmp_Row(2) = Tmp_Row(0);
    Tmp_Row(12) = Tmp_Row(0);
    Tmp_Row(14) = Tmp_Row(0);
    Tmp_Row(4) = -63.0 * zpx;
    Tmp_Row(6) = Tmp_Row(4);
    Tmp_Row(8) = Tmp_Row(4);
    Tmp_Row(10) = Tmp_Row(4);
    Tmp_Row(1) = -6.0 * zpx;
    Tmp_Row(13) = Tmp_Row(1);
    Tmp_Row(5) = -378.0 * zpx;
    Tmp_Row(9) = Tmp_Row(5);
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskZ.row(3).cols(0, 14) = Tmp_Row.cols(0, 14);
    EpsilonMaskZ.row(5).cols(1, 15) = Tmp_Row.cols(0, 14);

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = zpx;
    Tmp_Row(3) = Tmp_Row(0);
    Tmp_Row(12) = Tmp_Row(0);
    Tmp_Row(15) = Tmp_Row(0);
    Tmp_Row(4) = 63.0 * zpx;
    Tmp_Row(7) = Tmp_Row(4);
    Tmp_Row(8) = Tmp_Row(4);
    Tmp_Row(11) = Tmp_Row(4);
    Tmp_Row(1) = 7.0 * zpx;
    Tmp_Row(2) = Tmp_Row(1);
    Tmp_Row(13) = Tmp_Row(1);
    Tmp_Row(14) = Tmp_Row(1);
    Tmp_Row(5) = 441.0 * zpx;
    Tmp_Row(6) = Tmp_Row(5);
    Tmp_Row(9) = Tmp_Row(5);
    Tmp_Row(10) = Tmp_Row(5);
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskZ.row(4) = Tmp_Row;

    EpsilonMaskZ.save(hdf5_name("../Matrices.h5", "/EpsilonMaskMatrixZ", hdf5_opts::append));

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = -3.0 * xpz;
    Tmp_Row(2) = Tmp_Row(0);
    Tmp_Row(8) = Tmp_Row(0);
    Tmp_Row(10) = Tmp_Row(0);
    Tmp_Row(1) = -26.0 * xpz;
    Tmp_Row(9) = Tmp_Row(1);
    Tmp_Row(4) = -18.0 * xpz;
    Tmp_Row(6) = Tmp_Row(4);
    Tmp_Row(5) = -156.0 * xpz;
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskX.row(0).cols(0, 10) = Tmp_Row.cols(0, 10);
    EpsilonMaskX.row(2).cols(1, 11) = Tmp_Row.cols(0, 10);
    EpsilonMaskX.row(6).cols(4, 14) = Tmp_Row.cols(0, 10);
    EpsilonMaskX.row(8).cols(5, 15) = Tmp_Row.cols(0, 10);

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = -1.0 * xpz;
    Tmp_Row(3) = Tmp_Row(0);
    Tmp_Row(8) = Tmp_Row(0);
    Tmp_Row(11) = Tmp_Row(0);
    Tmp_Row(1) = -63.0 * xpz;
    Tmp_Row(2) = Tmp_Row(1);
    Tmp_Row(9) = Tmp_Row(1);
    Tmp_Row(10) = Tmp_Row(1);
    Tmp_Row(4) = -6.0 * xpz;
    Tmp_Row(7) = Tmp_Row(4);
    Tmp_Row(5) = -378.0 * xpz;
    Tmp_Row(6) = Tmp_Row(5);
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskX.row(1).cols(0, 11) = Tmp_Row.cols(0, 11);
    EpsilonMaskX.row(7).cols(4, 15) = Tmp_Row.cols(0, 11);

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = 3.0 * xpz;
    Tmp_Row(2) = Tmp_Row(0);
    Tmp_Row(12) = Tmp_Row(0);
    Tmp_Row(14) = Tmp_Row(0);
    Tmp_Row(4) = 21.0 * xpz;
    Tmp_Row(6) = Tmp_Row(4);
    Tmp_Row(8) = Tmp_Row(4);
    Tmp_Row(10) = Tmp_Row(4);
    Tmp_Row(1) = 26.0 * xpz;
    Tmp_Row(13) = Tmp_Row(1);
    Tmp_Row(5) = 182.0 * xpz;
    Tmp_Row(9) = Tmp_Row(5);
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskX.row(3).cols(0, 14) = Tmp_Row.cols(0, 14);
    EpsilonMaskX.row(5).cols(1, 15) = Tmp_Row.cols(0, 14);

    Tmp_Row.zeros(1, 16);
    Tmp_Row(0) = xpz;
    Tmp_Row(3) = Tmp_Row(0);
    Tmp_Row(12) = Tmp_Row(0);
    Tmp_Row(15) = Tmp_Row(0);
    Tmp_Row(4) = 7.0 * xpz;
    Tmp_Row(7) = Tmp_Row(4);
    Tmp_Row(8) = Tmp_Row(4);
    Tmp_Row(11) = Tmp_Row(4);
    Tmp_Row(1) = 63.0 * xpz;
    Tmp_Row(2) = Tmp_Row(1);
    Tmp_Row(13) = Tmp_Row(1);
    Tmp_Row(14) = Tmp_Row(1);
    Tmp_Row(5) = 441.0 * xpz;
    Tmp_Row(6) = Tmp_Row(5);
    Tmp_Row(9) = Tmp_Row(5);
    Tmp_Row(10) = Tmp_Row(5);
    Tmp_Row = Tmp_Row / 1536.0;
    EpsilonMaskX.row(4) = Tmp_Row;

    EpsilonMaskX.save(hdf5_name("../Matrices.h5", "/EpsilonMaskMatrixX", hdf5_opts::append));
    //-----------------------------------------------------------------------------------------------------------------------

    Mat<double> M_MaskMatrix;
    M_MaskMatrix.zeros(1, 9);
    M_MaskMatrix(0) = 1.0 / 36.0;
    M_MaskMatrix(2) = M_MaskMatrix(0);
    M_MaskMatrix(6) = M_MaskMatrix(0);
    M_MaskMatrix(8) = M_MaskMatrix(0);
    M_MaskMatrix(1) = 1.0 / 9.0;
    M_MaskMatrix(3) = M_MaskMatrix(1);
    M_MaskMatrix(5) = M_MaskMatrix(1);
    M_MaskMatrix(7) = M_MaskMatrix(1);
    M_MaskMatrix(4) = 4.0 / 9.0;
    M_MaskMatrix *= dX * dZ * 4.0 * pi * pi;
    M_MaskMatrix.save(hdf5_name("../Matrices.h5", "/M_MaskMatrix", hdf5_opts::append));

    //-----------------------------------------------------------------------------------------------------------------------
    Mat<double> TubeMaskZ,
        TubeMaskX;
    TubeMaskX.zeros(1, 3);
    TubeMaskX(0) = dX / 6.0;
    TubeMaskX(1) = 2.0 * dX / 3.0;
    TubeMaskX(2) = dX / 6.0;
    TubeMaskZ.zeros(Nz, 3);

    double Z = Zmin, D;

    for (int i = 0; i < Nz; i++)
    {
        Z += dZ;
        double atam = atan((dZ - Z) / L);
        double ata0 = atan(Z / L);
        double atap = atan((dZ + Z) / L);
        double logm = log((dZ - Z) * (dZ - Z) + L * L);
        double log0 = log(Z * Z + L * L);
        double logp = log((dZ + Z) * (dZ + Z) + L * L);
        /*         TubeMaskZ(i,0) = 0.0; (dZ - L*(atam+ata0))*2.0/dZ/dZ - 
                         Z*(log0-logm)/dZ/dZ;
        TubeMaskZ(i,1) = 0.0;(2.0*log0-logp-logm)/dZ +
                         (2.0*dZ - L*(atap+atam) )*2.0/dZ/dZ -
                         (logp+logm)*Z/dZ/dZ;   
        TubeMaskZ(i,2) = 0.0;(dZ - L*(atap-ata0))*2.0/dZ/dZ-
                         Z*(logp-log0)/dZ/dZ;    */
        D = (L * L + Z * Z);
        TubeMaskZ(i, 0) = -Z / D + (L * L - Z * Z) / D / D * dZ / 3.0 + (3.0 * L * L * Z - Z * Z * Z) / D / D / D * dZ * dZ / 6.0;
        TubeMaskZ(i, 1) = -2.0 * (L * L - Z * Z) / D / D * dZ;
        TubeMaskZ(i, 2) = +Z / D + (L * L - Z * Z) / D / D * dZ / 3.0 - (3.0 * L * L * Z - Z * Z * Z) / D / D / D * dZ * dZ / 6.0;

        //cout << TubeMaskZ(i,1) << endl;
    }
    TubeMaskX.save(hdf5_name("../Matrices.h5", "/Tube_XMaskMatrix", hdf5_opts::append));
    TubeMaskZ.save(hdf5_name("../Matrices.h5", "/Tube_ZMaskMatrix", hdf5_opts::append));

    //-----------------------------------------------------------------------------------------------------------------------
    Mat<int> NeibourN;
    NeibourN.zeros(Nmax + Nx, 9);

    for (int xi = 0; xi < Nx; xi++)
        for (int zi = 0; zi < Nz + 1; zi++)
        {
            ind = zi * Nx + xi;
            NeibourN(ind, 0) = (zi == 0) ? -1 : (xi % Nx == 0) ? ind - 1 : ind - Nx - 1;
            NeibourN(ind, 1) = (zi == 0) ? -1 : ind - Nx;
            NeibourN(ind, 2) = (zi == 0) ? -1 : (xi % Nx == Nx - 1) ? ind - 2 * Nx + 1 : ind - Nx + 1;
            NeibourN(ind, 3) = (xi % Nx == 0) ? ind + Nx - 1 : ind - 1;
            NeibourN(ind, 4) = ind;
            NeibourN(ind, 5) = (xi % Nx == Nx - 1) ? ind - Nx + 1 : ind + 1;
            NeibourN(ind, 6) = (zi == Nz) ? -1 : (xi % Nx == 0) ? ind + 2 * Nx - 1 : ind + Nx - 1;
            NeibourN(ind, 7) = (zi == Nz) ? -1 : ind + Nx;
            NeibourN(ind, 8) = (zi == Nz) ? -1 : (xi % Nx == Nx - 1) ? ind + 1 : ind + Nx + 1;
        }
    NeibourN.save(hdf5_name("../Matrices.h5", "/NeighboringConcElements", hdf5_opts::append));

    //-----------------------------------------------------------------------------------------------------------------------
    Mat<double> DiffMask;
    DiffMask.zeros(1, 9);
    DiffMask(4) = (zpx + xpz) * 4.0 / 3.0;
    DiffMask(1) = (zpx - 2.0 * xpz) / 3.0;
    DiffMask(7) = DiffMask(1);
    DiffMask(3) = (-2.0 * zpx + xpz) / 3.0;
    DiffMask(5) = DiffMask(3);
    DiffMask(0) = -(zpx + xpz) / 6.0;
    DiffMask(2) = DiffMask(0);
    DiffMask(6) = DiffMask(0);
    DiffMask(8) = DiffMask(0);
    DiffMask.save(hdf5_name("../Matrices.h5", "/Diff_MaskMatrix", hdf5_opts::append));
}