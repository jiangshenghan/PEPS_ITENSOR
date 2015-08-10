
#ifndef _PREDEF_H_
#define _PREDEF_H_

#include <iostream>
#include <complex>
#include <array>
#include <vector>
#include <algorithm>

using namespace std;
using namespace boost;
using namespace itensor;

enum Boundary {Open,Cylinder,Torus}


typedef std::complex<double> Complex;
//coordinate for both site and links, expressed as (x,y,s/b)
typedef array<int,3> Coordinate;

#endif
