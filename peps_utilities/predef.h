
#ifndef _PREDEF_H_
#define _PREDEF_H_

#include <iostream>
#include <sstream>
#include <complex>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <initializer_list>
#include <ctime>
#include <chrono>
#include <random>
#include <functional>
#include <core.h>
#include <wignerSymbols.h>
#include <armadillo>

//#define NDEBUG
#include <cassert>

using namespace itensor;
using namespace WignerSymbols;
//using namespace arma;

using std::cout;
using std::endl;

using Complex=std::complex<double>;
//coordinate for both sites and bonds, expressed as (x,y,s/b)
using Coordinate=std::array<int,3>;

static const double EPSILON=1E-10;
//static const Complex I()

#endif
