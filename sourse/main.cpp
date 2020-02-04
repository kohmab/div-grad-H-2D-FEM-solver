#include <iostream>
#include <armadillo>
#include <omp.h>

#include "./InitParams.h"
#include "./Absorber.h"
#include "./Source.h"
#include "./PlasProfile.h"
#include "./MyPhys.h"
#include "./IncBeam.h"
#include "./SetConfig.h"
#include "./EffectiveMedia.h"

using namespace arma;
using arma::span;
using arma::zeros;
using std::cin;
using std::complex;
using std::cout;

int main()

{
    ifstream ConfigFile;
    string ConfigString, SaveFileName = "";
    Mat<double> OptionsMat;
    OptionsMat.zeros(1, 23);

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

    int nthreads = 4;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, nthreads);
    omp_set_num_threads(nthreads);

    int PlotStep;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, PlotStep);

    int Nz = 2;
    int Nx = 2;

    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, Nz);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, Nx);
    int Nmax = Nx * Nz;
    SaveFileName += "Nx" + to_string(Nx) + "_Nz" + to_string(Nz);
    OptionsMat(0) = double(Nx);
    OptionsMat(1) = double(Nz);

    double dT = 1.0e-3, dTmin, dTmax, TmaxFs = 3.0, TFs = 0.0, T, Tmax, MaxNChange;
    Mat<double> Tm;
    Tm.zeros(1, 1);
    Tm[1] = T;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, dTmin);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, dTmax);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, MaxNChange);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, TFs);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, TmaxFs);
    dT = dTmax;
    T = TFs * 1.0e-15 * Field.Frequency / 2.0 / pi;
    Tmax = TmaxFs * 1.0e-15 * Field.Frequency / 2.0 / pi;
    SaveFileName += "_Tb" + d_to_string(TFs) + "_Te" + d_to_string(TmaxFs);
    OptionsMat(2) = dTmin;
    OptionsMat(3) = dTmax;

    double ZminLam, ZmaxLam;
    double XminLam, XmaxLam;
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
    SaveFileName += "_Z" + d_to_string(ZminLam) + "to" + d_to_string(ZmaxLam);
    SaveFileName += "_X" + d_to_string(XminLam) + "to" + d_to_string(XmaxLam);
    OptionsMat(5) = ZminLam;
    OptionsMat(6) = ZmaxLam;
    OptionsMat(7) = XminLam;
    OptionsMat(8) = XmaxLam;

    double AbWiLam, ReA, ImA, ReB, ImB;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, AbWiLam);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ReA);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ImA);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ReB);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ImB);
    SaveFileName += "_La" + d_to_string(AbWiLam);
    OptionsMat(9) = AbWiLam;

    double kLf;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, kLf);
    SaveFileName += "_kLf" + d_to_string(kLf);
    OptionsMat(10) = kLf;
    kLf = kLf / 2.0 / pi;

    int SourseType = 1;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, SourseType);
    double SourseDuFs, SourseShFs;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, SourseDuFs);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, SourseShFs);
    SaveFileName += "_tau" + d_to_string(SourseDuFs) + "_ts" + d_to_string(SourseShFs) + "_Xvar" + to_string(SourseType - 1);
    OptionsMat(11) = SourseDuFs * 1.0e-15 * Field.Frequency / 2.0 / pi;
    OptionsMat(12) = SourseShFs * 1.0e-15 * Field.Frequency / 2.0 / pi;

    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, Field.MaximumIntensity);
    SaveFileName += "_e12I" + d_to_string(Field.MaximumIntensity);
    Field.MaximumIntensity *= 1.0e12;
    OptionsMat(13) = Field.MaximumIntensity;

    double ImDe;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, ImDe);
    complex<double> De(0.0, ImDe);
    SaveFileName += "_nu" + d_to_string(ImDe);
    OptionsMat(14) = ImDe;

    double DiffusionCoof; // [cm^2/s]
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, DiffusionCoof);
    double a0; // [nm]
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, a0);
    double Diffn = DiffusionCoof / a0 / a0 * 1.0e14 / Field.Frequency * 2.0 * pi;
    SaveFileName += "_D" + d_to_string(DiffusionCoof);
    SaveFileName += "_a" + d_to_string(a0);

    OptionsMat(15) = DiffusionCoof;

    double f, apb;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, f);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, apb);
    OptionsMat(20) = f;
    OptionsMat(21) = apb;
    OptionsMat(22) = a0;
    double Nc = apb * apb * f;

    SaveFileName += "_f" + d_to_string(f) + "_ab" + d_to_string(apb);

    string SilicaInside;
    int PertNx, PertNz;
    double PertAmplitude, StartHomConc;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, SilicaInside);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, PertNx);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, PertNz);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, PertAmplitude);
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, StartHomConc);
    OptionsMat(16) = double(PertNx);
    OptionsMat(17) = double(PertNz);
    OptionsMat(18) = PertAmplitude;
    OptionsMat(19) = StartHomConc;

    double SigmaAIns = 0.01, SigmaMIns = 0.01, RecFr = 0.01, Nat0, AbsE2Th;
    Nat0 = Silica.AtomConcentration / Field.CriticalConcentration / Silica.Permittivity;
    RecFr = 2.0 * pi / Silica.RecombinationTime / Field.Frequency;
    double AlphaFsIns = 8.0 * pi / 3.0 * Electron.Charge * Electron.Charge /
                        sqrt(Silica.Permittivity) / Electron.Mass / LightSpeed / Silica.BandGapInside / Field.Frequency * 1.0e7; /* [cm^2/W] */
    SigmaMIns = 2.0 * pi * Silica.SigmaInside * pow(Field.MaximumIntensity, Silica.PhotonsNumInside) /
                Field.Frequency * Nat0;
    SigmaAIns = 2.0 * pi * AlphaFsIns * Field.MaximumIntensity / Field.Frequency / Silica.Permittivity * ImDe / (1.0 + ImDe * ImDe);
    double CollFrAlpha = Field.CriticalConcentration * Silica.BandGapInside * LightSpeed * sqrt(Silica.Permittivity) /
                         Field.MaximumIntensity / 1.0e7;

    double SigmaAOut = 0.01, SigmaMOut = 0.01;
    double AlphaFsOut = 8.0 * pi / 3.0 * Electron.Charge * Electron.Charge /
                        sqrt(Silica.Permittivity) / Electron.Mass / LightSpeed / Silica.BandGapOutside / Field.Frequency * 1.0e7; /* [cm^2/W] */
    SigmaMOut = 2.0 * pi * Silica.SigmaOutside * pow(Field.MaximumIntensity, Silica.PhotonsNumOutside) /
                Field.Frequency * Nat0;
    SigmaAOut = 2.0 * pi * AlphaFsOut * Field.MaximumIntensity / Field.Frequency / Silica.Permittivity * ImDe / (1.0 + ImDe * ImDe);

    SaveFileName = SilicaInside + "_" + SaveFileName;
    SaveFileName += "_Pert_Nx" + to_string(PertNx) + "_Nz" + to_string(PertNz) +
                    "_A" + d_to_string(PertAmplitude) + "_Ns" + d_to_string(OptionsMat(19));
    int Power = Silica.PhotonsNumInside;

    if (SilicaInside == "yes")
    {
        Power = Silica.PhotonsNumOutside;
        SigmaAIns = SigmaAOut;
        SigmaMIns = SigmaMOut;
        CollFrAlpha = Field.CriticalConcentration * Silica.BandGapOutside * LightSpeed * sqrt(Silica.Permittivity) /
                      Field.MaximumIntensity / 1.0e7;
    }

    AbsE2Th = RecFr / SigmaAIns;

    OptionsMat(4) = AbsE2Th;

    double FieldCalibrationConst;
    getline(ConfigFile, ConfigString);
    GetParam(ConfigString, FieldCalibrationConst);

    ConfigFile.close();
    SaveFileName = "../Data/" + SaveFileName + ".h5";
    cout << "Data will be saved in " << SaveFileName << endl;
    OptionsMat.save(hdf5_name(SaveFileName, "/Options"));

    cx_colvec Zh = linspace<cx_colvec>(Zmin, Zmax, Nz),
              Zn = linspace<cx_colvec>(Zmin - dZ / 2.0, Zmax + dZ / 2.0, Nz + 1),
              Z_plot = Zh * optp, Nz_plot = linspace<cx_colvec>(0.0, double(Nz - 1), Nz);
    cx_rowvec X = linspace<cx_rowvec>(Xmin, Xmax, Nx),
              X_plot = X * optp, Nx_plot = linspace<cx_rowvec>(0.0, double(Nx - 1), Nx);

    cx_mat Z2D;
    Z2D.zeros(Nz, Nx);
    Z2D.each_col() += Zh;
    cx_mat X2D;
    X2D.zeros(Nz, Nx);
    X2D.each_row() += X;

    cx_vec plt_vec, TubeFunc;

    MyAbsorber Absorber(Zn, Zh, AbWiLam, ReA + I * ImA, ReB + I * ImB);

    MySource Source(SourseType, Zh, X, kLf, Absorber.Width,
                    SourseDuFs * 1.0e-15 * Field.Frequency / 2.0 / pi, SourseShFs * 1.0e-15 * Field.Frequency / 2.0 / pi, FieldCalibrationConst);

    /*     cout << Me.eps_x(cx_double(0.0)) << " " << Me.eps_z(cx_double(0.0)) << endl;
    return 1; */
    GenRayTubeProfile(TubeFunc, Zh, kLf, Absorber.Width, Source.ZWidth, Nz);
    GenRayTubeFunc(TubeFunc);

    // CalcSourceAmplitude(Source, Absorber, TubeFunc, dZ);
    Source.Amplitude = Source.Amplitude / sqrt(1.0 + (Zmin + Source.ZShift + Source.ZWidth) * (Zmin + Source.ZShift + Source.ZWidth) / kLf / kLf);
    Source.Shape *= Source.Amplitude;

    // GenAtomDistribution(Nat, Z, Nat0, Absorber.Width, Source.ZWidth);

    Mat<cx_double> H2D(Nz, Nx, fill::zeros),
        Hy0(Nmax, 1, fill::zeros), Hy1(Nmax, 1, fill::zeros),
        Ex0(Nmax + Nx, 1, fill::zeros), Ex1(Nmax + Nx, 1, fill::zeros),
        Ez0(Nmax, 1, fill::zeros), Ez1(Nmax, 1, fill::zeros),
        Exa(Nmax + Nx, 1, fill::zeros), Eza(Nmax, 1, fill::zeros),
        EpsX(Nmax + Nx, 1, fill::zeros), EpsZ(Nmax + Nx, 1, fill::zeros), Eps2D;

    Mat<double> AbsE2(Nmax + Nx, 1, fill::zeros),
        N0(Nmax + Nx, 1, fill::zeros), N05(Nmax + Nx, 1, fill::zeros), N1(Nmax + Nx, 1, fill::zeros),
        Nat(Nmax + Nx, 1, fill::zeros), EpsM(Nmax + Nx, 1, fill::zeros),
        Heat(Nmax + Nx, 1, fill::zeros), CollFr(Nmax + Nx, 1, fill::zeros);

    int plt = 0,
        sav = 0;
    double ExecTime = omp_get_wtime();
    cout << "dX = " << dX << ", dZ = " << dZ << endl;

    // Load the prepeared matrices
    Mat<int> H_ind, Eps_ind, H_neib, Eps_neib, N_neib;
    Mat<double> EpsX_mask, EpsZ_mask, M_mask, T_zmask, T_xmask, D_mask;

    H_ind.load(hdf5_name("../Matrices.h5", "ElementsIndexes"));
    Eps_ind.load(hdf5_name("../Matrices.h5", "PermittivityIndexes"));
    H_neib.load(hdf5_name("../Matrices.h5", "NeighboringElements"));
    Eps_neib.load(hdf5_name("../Matrices.h5", "PermittivityAround"));
    N_neib.load(hdf5_name("../Matrices.h5", "NeighboringConcElements"));

    EpsX_mask.load(hdf5_name("../Matrices.h5", "EpsilonMaskMatrixX"));
    EpsZ_mask.load(hdf5_name("../Matrices.h5", "EpsilonMaskMatrixZ"));
    M_mask.load(hdf5_name("../Matrices.h5", "M_MaskMatrix"));

    T_zmask.load(hdf5_name("../Matrices.h5", "Tube_ZMaskMatrix"));
    T_xmask.load(hdf5_name("../Matrices.h5", "Tube_XMaskMatrix"));

    D_mask.load(hdf5_name("../Matrices.h5", "Diff_MaskMatrix"));

    // Configuring the perturbation;
    MyPerturbation Pert(Nx, Nz + 1, PertNx, PertNz, PertAmplitude, X, Zn, kLf);

    double tm = omp_get_wtime();

    bool IonizationGoes = false;

    for (uword i = 0; i < Nmax + Nx; i++)
    {
        N0(i) = (Pert.ZShape(Eps_ind(i, 0)) + Pert.Shape(Eps_ind(i, 0), Eps_ind(i, 1))) * StartHomConc;
    }

    N1 = N0;

    MyMedia Me(Nmax + Nx, apb, f, Nc, De, Diffn, SigmaAIns, N0);
    for (uword i = 0; i < Nmax + Nx; i++)
    {
        EpsX(i) = Me.eps_x(i) - Absorber.ShapeN(Eps_ind(i, 0));
        EpsZ(i) = Me.eps_z(i) - Absorber.ShapeN(Eps_ind(i, 0));
    }

    EpsX.save(hdf5_name(SaveFileName, "/EpsX_at_step_0", hdf5_opts::append));
    EpsZ.save(hdf5_name(SaveFileName, "/EpsZ_at_step_0", hdf5_opts::append));
    Pert.Shape.save(hdf5_name(SaveFileName, "/Pert", hdf5_opts::append));

    // Generating Matrices;
    Mat<cx_double> S = zeros<Mat<cx_double>>(Nmax, 1);

    umat LocationsH, LocationsN;
    Col<cx_double> Mv, Rv, Uv;
    LocationsH.zeros(2, ((Nz - 2) * 9 + 12) * Nx);
    LocationsN.zeros(2, (Nz * 9 + 12) * Nx);
    Mv.zeros(((Nz - 2) * 9 + 12) * Nx);
    Rv.zeros(((Nz - 2) * 9 + 12) * Nx);
    Uv.zeros(((Nz - 2) * 9 + 12) * Nx);

    uword spind = 0;
    //#pragma omp parallel for
    for (uword i = 0; i < Nmax; i++)
    {
        Mat<cx_double> EpsX_vec = zeros<Mat<cx_double>>(16, 1),
                       EpsZ_vec = zeros<Mat<cx_double>>(16, 1),
                       Eps_H_Mask = zeros<Mat<cx_double>>(9, 1);
        for (uword j = 0; j < 16; j++)
        {
            if (Eps_neib(i, j) >= 0)
            {
                EpsX_vec(j) = 1.0 / EpsX(Eps_neib(i, j));
                EpsZ_vec(j) = 1.0 / EpsZ(Eps_neib(i, j));
            }
        }
        Eps_H_Mask = EpsX_mask * EpsX_vec + EpsZ_mask * EpsZ_vec;

        for (uword j = 0; j < 9; j++)
        {
            if (H_neib(i, j) >= 0)
            {
                LocationsH(0, spind) = i;
                LocationsH(1, spind) = H_neib(i, j);
                Mv(spind) = M_mask(j) *
                            (1.0 - Absorber.ShapeH(H_ind(H_neib(i, j), 0)));
                Rv(spind) = Eps_H_Mask(j);
                Uv(spind) = T_xmask(j % 3) * T_zmask(H_ind(i, 0), j / 3) * (1.0 - Absorber.ShapeH(H_ind(H_neib(i, j), 0)));
                S(i) += M_mask(j) * Source.Shape(H_ind(H_neib(i, j), 0), H_ind(H_neib(i, j), 1)) * I * pi;
                spind++;
            }
        }
    }

    SpMat<cx_double> M(LocationsH, Mv, Nmax, Nmax, true, true);
    SpMat<cx_double> R(LocationsH, Rv, Nmax, Nmax, true, true);
    SpMat<cx_double> U(LocationsH, Uv, Nmax, Nmax, true, true);

    cout << "Matrix initialization needs " << omp_get_wtime() - tm << " seconds." << endl;

    SpMat<cx_double> Denum = (1.0 + pi * I * dT / 2.0) * M - pi * I * dT / 2.0 * (R - U),
                     Num = (1.0 - pi * I * dT / 2.0) * M + pi * I * dT / 2.0 * (R - U);
    /* ------------------------------------------ */
    //Source.Shape.zeros(Nz);
    //Hy00 = exp(-Z2D%Z2D + -X2D%X2D*100.0);
    /* ------------------------------------------ */
    Tm[0] = T;
    Tm.save(hdf5_name(SaveFileName, "/Time_at_step_" + to_string(0), hdf5_opts::append));
    Hy1.save(hdf5_name(SaveFileName, "/H_at_step_" + to_string(0), hdf5_opts::append));
    Ex1.save(hdf5_name(SaveFileName, "/Ex_at_step_" + to_string(0), hdf5_opts::append));
    Ez1.save(hdf5_name(SaveFileName, "/Ez_at_step_" + to_string(0), hdf5_opts::append));
    AbsE2.save(hdf5_name(SaveFileName, "/|E|2_at_step_" + to_string(0), hdf5_opts::append));
    N0.save(hdf5_name(SaveFileName, "/N_at_step_" + to_string(0), hdf5_opts::append));
    Heat.save(hdf5_name(SaveFileName, "/Heat_at_step_" + to_string(0), hdf5_opts::append));
    Me.f.save(hdf5_name(SaveFileName, "/Fill_at_step_" + to_string(0), hdf5_opts::append));
    Me.apb.save(hdf5_name(SaveFileName, "/SzRatio_at_step_" + to_string(0), hdf5_opts::append));

    while (T < Tmax)
    {
        // Find H_y;
        Hy1 = spsolve(Denum, Num * Hy0 + S * Source.TimeFunctionIntegral(T, dT));
        // Find E_x;
        for (int i = 0; i < Nx; i++)
        {
            Ex1(i) = ((1.0 - pi * I * dT) * Ex0[i] -
                      dT / (EpsX[Eps_neib(i, 5)] + EpsX[Eps_neib(i, 6)]) *
                          (Hy1[i] + Hy0[i]) / dZ) /
                     (1.0 + pi * I * dT);
        }
#pragma omp parallel for
        for (int i = Nx; i < Nmax; i++)
        {
            Ex1(i) = ((1.0 - pi * I * dT) * Ex0(i) -
                      dT / (EpsX[Eps_neib(i, 5)] + EpsX[Eps_neib(i, 6)]) *
                          ((Hy1[i] + Hy0[i]) - (Hy1[i - Nx] + Hy0[i - Nx])) / dZ) /
                     (1.0 + pi * I * dT);
        }
        for (int i = Nmax; i < Nmax + Nx; i++)
        {
            Ex1(i) = ((1.0 - pi * I * dT) * Ex0[i] -
                      dT / (EpsX[Eps_neib(i - Nx, 10)] + EpsX[Eps_neib(i - Nx, 9)]) *
                          (-(Hy1[i - Nx] + Hy0[i - Nx])) / dZ) /
                     (1.0 + pi * I * dT);
        }
// Find E_z;
#pragma omp parallel for
        for (int i = 0; i < Nmax; i++)
        {
            Ez1(i) = ((1.0 - pi * I * dT) * Ez0[i] +
                      dT / (EpsZ[Eps_neib(i, 6)] + EpsZ[Eps_neib(i, 10)]) *
                          ((Hy1[H_neib(i, 5)] + Hy0[H_neib(i, 5)]) - (Hy1[i] + Hy0[i])) / dX) /
                     (1.0 + pi * I * dT);
        }
        // Find |E|^2;

#pragma omp for
        for (int i = 0; i < Nmax; i++)
        {
            Exa[i] = (Ex1[i] + Ex1[H_neib(i, 5)]) / 2.0;
        }
        for (int i = Nmax; i < Nmax + Nx; i++)
        {
            Exa[i] = (Ex1[Eps_neib(i - Nx, 10)] + Ex1[Eps_neib(i - Nx, 11)]) / 2.0;
        }
        for (int i = 0; i < Nx; i++)
        {
            Eza[i] = Ez1[i] / 2.0;
        }
#pragma omp for
        for (int i = Nx; i < Nmax + Nx; i++)
        {
            Eza[i] = (Ez1[i] + Eza[i - Nx]) / 2.0;
        }
#pragma omp parallel for
        for (int i = 0; i < Nmax + Nx; i++)
        {
            AbsE2[i] = real(Eza[i] * Me.al_z[i] * conj(Eza[i] * Me.al_z[i]) +
                            Exa[i] * Me.al_x[i] * conj(Exa[i] * Me.al_x[i]));
        }
        if (!IonizationGoes)
        {
            if (AbsE2.max() > AbsE2Th)
                IonizationGoes = true;
        }
        // Find N;
        if (IonizationGoes)
        {
#pragma omp parallel for
            for (uword i = 0; i < Nmax + Nx; i++)
            {
                N1[i] = ((1.0 + dT * ((Nat0 - N0[i]) / Nat0 * SigmaAIns * AbsE2[i] - RecFr)) * N0[i] + dT * (Nat0 - N0[i]) / Nat0 * SigmaMIns * pow(AbsE2[i], Power)) /
                        (1.0 - dT * ((Nat0 - N0[i]) / Nat0 * SigmaAIns * AbsE2[i] - RecFr));
                N1[i] = (N1[i] > Nat0 || N1[i] < 0.0) ? Nat0 - MaxNChange : N1[i];
            }
            Me.update(N0, N1, Exa, Eza, dT);
#pragma omp parallel for
            for (uword i = 0; i < Nmax + Nx; i++)
            {
                EpsX(i) = Me.eps_x(i) - Absorber.ShapeN(Eps_ind(i, 0));
                EpsZ(i) = Me.eps_z(i) - Absorber.ShapeN(Eps_ind(i, 0));
            }
            // Gen new R-matrix
            spind = 0;
            for (uword i = 0; i < Nx; i++)
            {
                Mat<cx_double> EpsX_vec = zeros<Mat<cx_double>>(16, 1),
                               EpsZ_vec = zeros<Mat<cx_double>>(16, 1),
                               Eps_H_Mask = zeros<Mat<cx_double>>(9, 1);
                for (uword j = 0; j < 16; j++)
                {
                    if (Eps_neib(i, j) >= 0)
                    {
                        EpsX_vec(j) = 1.0 / EpsX(Eps_neib(i, j));
                        EpsZ_vec(j) = 1.0 / EpsZ(Eps_neib(i, j));
                    }
                }
                Eps_H_Mask = EpsX_mask * EpsX_vec + EpsZ_mask * EpsZ_vec;

                for (uword j = 0; j < 9; j++)
                {
                    if (H_neib(i, j) >= 0)
                    {
                        Rv(spind) = Eps_H_Mask(j);
                        spind++;
                    }
                }
            }
#pragma omp parallel for
            for (uword i = Nx; i < Nmax - Nx; i++)
            {
                Mat<cx_double> EpsX_vec = zeros<Mat<cx_double>>(16, 1),
                               EpsZ_vec = zeros<Mat<cx_double>>(16, 1),
                               Eps_H_Mask = zeros<Mat<cx_double>>(9, 1);
                for (uword j = 0; j < 16; j++)
                {
                    EpsX_vec(j) = 1.0 / EpsX(Eps_neib(i, j));
                    EpsZ_vec(j) = 1.0 / EpsZ(Eps_neib(i, j));
                }
                Eps_H_Mask = EpsX_mask * EpsX_vec + EpsZ_mask * EpsZ_vec;

                for (uword j = 0; j < 9; j++)
                    Rv(6 * Nx + j + 9 * (i - Nx)) = Eps_H_Mask(j);
            }
            spind = 6 * Nx + 9 * Nx * (Nz - 2);
            for (uword i = Nmax - Nx; i < Nmax; i++)
            {
                Mat<cx_double> EpsX_vec = zeros<Mat<cx_double>>(16, 1),
                               EpsZ_vec = zeros<Mat<cx_double>>(16, 1),
                               Eps_H_Mask = zeros<Mat<cx_double>>(9, 1);
                for (uword j = 0; j < 16; j++)
                {
                    if (Eps_neib(i, j) >= 0)
                    {
                        EpsX_vec(j) = 1.0 / EpsX(Eps_neib(i, j));
                        EpsZ_vec(j) = 1.0 / EpsZ(Eps_neib(i, j));
                    }
                }
                Eps_H_Mask = EpsX_mask * EpsX_vec + EpsZ_mask * EpsZ_vec;

                for (uword j = 0; j < 9; j++)
                {
                    if (H_neib(i, j) >= 0)
                    {
                        Rv(spind) = Eps_H_Mask(j);
                        spind++;
                    }
                }
            }
            SpMat<cx_double> Rnew(LocationsH, Rv, Nmax, Nmax, true, true);
            R = Rnew;
            Denum = (1.0 + pi * I * dT / 2.0) * M - pi * I * dT / 2.0 * (R - U);
            Num = (1.0 - pi * I * dT / 2.0) * M + pi * I * dT / 2.0 * (R - U);
        }
        // ----------
        if (abs((N1.max() - N0.max())) <= MaxNChange || dT == dTmin || AbsE2.max() < AbsE2Th)
        {
            if (dT < dTmax && AbsE2.max() < AbsE2Th && T > Source.TShift)
            {
                dT = dTmax;
                cout << "Time step is changed to dT = " << dT << endl;
            }
            else
            {
#pragma omp parallel for
                for (int i = 0; i < Nmax + Nx; i++)
                {
                    CollFr(i) =
                        (N1(i) > N0(i) && AbsE2(i) > 1.0e-3 && N1(i) > 1.0e-4)
                            ? imag(De) + (0.5 + CollFrAlpha / AbsE2[i]) * (N1(i) - N0(i)) / (N1(i) + N0(i)) / M_PI / dT
                            : imag(De);
                    Heat(i) = Heat(i) + Me.f[i] * 0.25 / M_PI * AbsE2(i) * CollFr(i) * (N1(i) + N0(i)) * dT;
                }

                Hy0 = Hy1;
                Ex0 = Ex1;
                Ez0 = Ez1;
                N0 = N1;
                T += dT;
                plt++;
            }
        }
        else
        {
            double dTtmp = dT * 0.75 * MaxNChange / abs(N1.max() - N0.max());
            dT = dTtmp > dTmin ? dTtmp : dTmin;
            cout << "Time step is changed to dT = " << dT << endl;
        }
        /* ------------------------------------------ */
        if (plt % PlotStep == 0)
        {

            cout << PlotStep << " steps need " << omp_get_wtime() - ExecTime << " seconds." << endl;
            cout << "Amplitude at T = " << T << " is " << Source.TimeFunction(T + dT / 2.0)
                 << "; dT = " << dT
                 << "; Tmax = " << Tmax << endl;

            plt = 0;

            sav++;
            Tm[0] = T;
            Tm.save(hdf5_name(SaveFileName, "/Time_at_step_" + to_string(sav), hdf5_opts::append));
            Hy1.save(hdf5_name(SaveFileName, "/H_at_step_" + to_string(sav), hdf5_opts::append));
            Ex1.save(hdf5_name(SaveFileName, "/Ex_at_step_" + to_string(sav), hdf5_opts::append));
            Ez1.save(hdf5_name(SaveFileName, "/Ez_at_step_" + to_string(sav), hdf5_opts::append));
            AbsE2.save(hdf5_name(SaveFileName, "/|E|2_at_step_" + to_string(sav), hdf5_opts::append));
            N1.save(hdf5_name(SaveFileName, "/N_at_step_" + to_string(sav), hdf5_opts::append));
            EpsX.save(hdf5_name(SaveFileName, "/EpsX_at_step_" + to_string(sav), hdf5_opts::append));
            EpsZ.save(hdf5_name(SaveFileName, "/EpsZ_at_step_" + to_string(sav), hdf5_opts::append));
            Heat.save(hdf5_name(SaveFileName, "/Heat_at_step_" + to_string(sav), hdf5_opts::append));
            Me.f.save(hdf5_name(SaveFileName, "/Fill_at_step_" + to_string(sav), hdf5_opts::append));
            Me.apb.save(hdf5_name(SaveFileName, "/SzRatio_at_step_" + to_string(sav), hdf5_opts::append));

            ExecTime = omp_get_wtime();
        }

        if (AbsE2.has_nan() || AbsE2.has_inf())
        {
            cout << "ERROR! Field is infinite or undefined.";
            return 1;
        }
    }
    return 0;
}