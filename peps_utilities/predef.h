
#ifndef _PREDEF_H_
#define _PREDEF_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <numeric>
#include <initializer_list>
#include <functional>
#include <ctime>
#include <chrono>
#include <random>
#include <functional>
#include <armadillo>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multimin.h>
#include "mpi.h"
#include "core.h"

//#define NDEBUG
#include <cassert>

using namespace itensor;
//using namespace arma;

using std::cout;
using std::endl;

using Complex=std::complex<double>;
//coordinate for both sites and bonds, expressed as (x,y,s/b)
using Coordinate=std::array<int,3>;

static const double EPSILON=1E-15;
//static const Complex I()

//random number generator between (-1,1)
static std::default_random_engine generator(std::time(0));
static std::uniform_real_distribution<double> distribution(-1.0,1.0);
static auto rand_gen=std::bind(distribution,generator);

#endif
