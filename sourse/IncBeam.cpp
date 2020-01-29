
#include "./IncBeam.h"

using namespace std;
using namespace arma;

void GenRayTubeProfile(cx_vec &S, cx_vec const &Z, double kLf, double AbsorberWidth, double SourceWidth, int Nz)
{
    S.zeros(Nz);
    double Zmin = real(Z(0)), Zmax = real(Z(Nz - 1)), L = Zmax - Zmin,
           buf = 2.0 * M_PI, fact = 1.0 / 5.0;
    int minind = floor(double(Nz) * (AbsorberWidth + SourceWidth) / L),
        maxind = Nz - floor(double(Nz) * AbsorberWidth / L),
        minind_buf = floor(double(Nz) * (AbsorberWidth + SourceWidth + buf) / L),
        maxind_buf = Nz - floor(double(Nz) * (AbsorberWidth + buf) / L);

    /*     for (int i = 0; i < minind; i++) 
        S(i) = (1.0 + (Z(minind)/kLf)*(Z(minind)/kLf)); */
    for (int i = 0; i < Nz; i++)
        S(i) = (1.0 + (Z(i) / kLf) * (Z(i) / kLf));
    /*     for (int i = maxind; i < Nz; i++) 
        S(i) = (1.0 + (Z(maxind)/kLf)*(Z(maxind)/kLf)); */

    /*     double a = real(S(minind))*fact + (1.0-fact)*real(S(minind_buf)), 
           b = real(S(minind_buf)), Lind = double(minind_buf - minind),
           alpha = real(S(minind_buf - 1) - S(minind_buf + 1)) / 2.0,
           A0 = a, A2 = 3.0*(b-a)/Lind/Lind+alpha/Lind , A3 = (2.0*(a - b) - Lind*alpha) / Lind / Lind / Lind;

    for (int i = minind; i < minind_buf + 1; i++) 
        S(i) = A0 + 
               A2 * double(i - minind) * double(i - minind) +    
               A3 * double(i - minind) * double(i - minind) * double(i - minind);

    a = real(S(maxind))*fact + (1.0-fact)*real(S(maxind_buf)),
    b = real(S(maxind_buf)); Lind = double(maxind - maxind_buf);
    alpha = real(S(maxind_buf + 1) - S(maxind_buf - 1)) / 2.0;
    A0 = a; A2 = 3.0*(b-a)/Lind/Lind+alpha/Lind; A3 = (2.0*(a - b) - Lind*alpha) / Lind / Lind / Lind;

    for (int i = maxind_buf; i < maxind + 1; i++) 
        S(i) = A0 + 
               A2 * double(maxind - i) * double(maxind - i) +    
               A3 * double(maxind - i) * double(maxind - i) * double(maxind - i);

    for (int i = 0; i < minind; i++) 
        S(i) = S(minind);
    for (int i = maxind; i< Nz; i++) 
        S(i) = S(maxind); */
}

void GenRayTubeFunc(cx_vec &S)
{
    /*   double Zmin = real(Z(0)), Zmax = real(Z(Nz-1)), L = Zmax - Zmin; 
    int minind = floor(double(Nz)*AbsorberWidth/L) , maxind = Nz - minind;     */
    cx_vec Tmp = S;
    S(0) = 0;
    for (int i = 1; i < S.n_elem - 1; i++)
        S(i) = (Tmp(i + 1) - Tmp(i - 1)) / 2.0 / Tmp(i);
    S(S.n_elem - 1) = 0;
}
