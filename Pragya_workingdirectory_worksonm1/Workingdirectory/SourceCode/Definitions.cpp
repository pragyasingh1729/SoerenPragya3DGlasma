#ifndef __DEFINITIONS_CPP__
#define __DEFINITIONS_CPP__

#include <ctime>
#include <cmath>
#include <math.h>
#include <cstring>
#include <cstdlib>
#include <ccomplex>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <vector>
// #include <x86intrin.h>

#ifdef __x86_64__
 #include <immintrin.h>
#else   
 #include "sse2neon.h"
#endif

#include <unistd.h>

//DEFINITIONS OF DATATYPES
typedef int INT;
typedef double DOUBLE;
typedef std::complex<DOUBLE> COMPLEX;

//DETERMINE MAXIMUM WORKING PRECISION
static const int MAX_DIGITS_PRECISION=std::numeric_limits<DOUBLE>::digits10;
static const int OUTPUT_PRECISION=MAX_DIGITS_PRECISION;

using std::sqrt;
using std::pow;
using std::abs;
using std::cos;
using std::sin;
using std::exp;
using std::atan2;

//BASIC CONSTANTS
//#ifndef ComplexI
#define ComplexI COMPLEX(0.0,1.0) // i
//#endif

//BASIC MACROS
#define SQR(x)      ((x)*(x))                        // x^2 
#define SQR_ABS(x)  (SQR(real(x)) + SQR(imag(x)))  // |x|^2
#define MOD(x,y)    ((x)>=(0)?((x)%(y)):(((x)+(y))%(y)))       // x mod y
#define DELTA(x,y)    (((x)==(y))?(1):(0))       // Kronecker delta
#define SIGN(x)  ((x)>0?(1):((x)<0?(-1):(0))) // SIGN FUNCTION

//COMPILER FLAGS FOR SUNc GAUGE GROUP
#define U1_FLAG  111
#define SU2_FLAG 222
#define SU3_FLAG 333

// COMPILER FLAGS FOR BOUNDARY CONDITIONS
#define OPEN_BOUNDARY_CONDITIONS_FLAG 444
#define PERIODIC_BOUNDARY_CONDITIONS_FLAG 555

// OPEN MP //
#include <omp.h>

// FFTW INCLUSION //
#include "fftw3.h"

#ifndef MY_FFTW_PLANNER_FLAG
#define MY_FFTW_PLANNER_FLAG FFTW_ESTIMATE
#endif

#endif
