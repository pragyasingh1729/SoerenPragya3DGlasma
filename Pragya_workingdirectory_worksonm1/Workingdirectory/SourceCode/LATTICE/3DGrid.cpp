#ifndef __LATTICE__CPP__
#define __LATTICE__CPP__

///////////////////////////////////////////////
//CONTAINS LATTICE DISCRETIZATION PARAMETERS //
///////////////////////////////////////////////

// CHOOSE PERIODIC OR OPEN BOUNDARY CONDTIONS

#define BOUNDARY_CONDITIONS_FLAG OPEN_BOUNDARY_CONDITIONS_FLAG

//#define BOUNDARY_CONDITIONS_FLAG PERIODIC_BOUNDARY_CONDITIONS_FLAG

namespace Lattice{
    
    //DIMENSION OF THE LATTICE
    static const int Dimension=3;

    //NUMBER OF LATTICE SITES //
    INT N[3]={128,128,2048};
    
    // LATTICE SPACINGS //0.01171 //NEXT TIME WHILE PUTTING VALUE FOR AZ FOR REALISTIC CASE RE CALCULATE IT?/
    DOUBLE a[3]={1.0,1.0,0.25};
    
    // SCALE //
    DOUBLE aScale=1.0;
    
    // VOLUME //
    INT Volume;
    
    // CUBIC SPACING //
    DOUBLE aCube;

}

#endif
