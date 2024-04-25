#ifndef __FOURIER_SPACE_CPP__
#define __FOURIER_SPACE_CPP__

// INCLUDE THREE DIMENSIONAL FAST FOURIER TRANSFORM //

#include "../FFT/FFT3D.cpp"
#include "../FFT/FFT2D.cpp"
#include "../FFT/FFT1D.cpp"
namespace FourierSpace{

    FFT3D *GaugeUpdate;

    // INITIALIZATION //
    void Init(){

	   GaugeUpdate=new FFT3D(SUNcGaugeLinks::U->N[0],SUNcGaugeLinks::U->N[1],SUNcGaugeLinks::U->N[2],SUNcAlgebra::NumberOfGenerators);
    }
}

#endif

