#pragma once
#ifndef SETCONFIG_H // include guard
#define SETCONFIG_H

#include <fstream>
#include <iostream>

using namespace std;

void GetParam(string &St, int &Param);
void GetParam(string &St, double &Param);
void GetParam(string &St, string &Param);
string d_to_string(double &St);

#endif
