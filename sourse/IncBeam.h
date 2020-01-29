#pragma once
#ifndef CALCEXTFIELD_H
#define CALCEXTFIELD_H

#include <armadillo>
//#include "./MyPhys.h"

using namespace arma;

void GenRayTubeProfile(cx_vec &S, cx_vec const &Z, double kLf, double AbsorberWidth, double SourceWidth, int Nz);
void GenRayTubeFunc(cx_vec &S);

#endif