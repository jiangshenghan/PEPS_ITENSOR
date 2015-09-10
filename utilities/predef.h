
#ifndef _PREDEF_H_
#define _PREDEF_H_

#include <iostream>
#include <sstream>
#include <complex>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <initializer_list>
#include <cassert>
#include "core.h"

using namespace itensor;

using std::cout;
using std::endl;

using Complex=std::complex<double>;
//coordinate for both sites and bonds, expressed as (x,y,s/b)
using Coordinate=std::array<int,3>;

#endif