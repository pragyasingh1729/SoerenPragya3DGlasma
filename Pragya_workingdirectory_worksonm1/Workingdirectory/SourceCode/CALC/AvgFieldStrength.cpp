#ifndef __AVG_FIELD_STRENGTH__CPP__
#define __AVG_FIELD_STRENGTH__CPP__

///////////////////////////////////////////////////
// COMPUTATION OF ADJOINT FOR PARALLEL TRANSPORT //
///////////////////////////////////////////////////

#define SET_ADJOINT_BUFFERS() \
SU_Nc_ADJOINT_FORMAT U0Adjx[SUNcAlgebra::NumberOfGenerators*SUNcAlgebra::NumberOfGenerators];\
SU_Nc_ADJOINT_FORMAT U1Adjy[SUNcAlgebra::NumberOfGenerators*SUNcAlgebra::NumberOfGenerators];\
SU_Nc_ADJOINT_FORMAT U2Adjz[SUNcAlgebra::NumberOfGenerators*SUNcAlgebra::NumberOfGenerators];\
\
SU_Nc_ADJOINT_FORMAT U0AdjxM[SUNcAlgebra::NumberOfGenerators*SUNcAlgebra::NumberOfGenerators];\
SU_Nc_ADJOINT_FORMAT U1AdjyM[SUNcAlgebra::NumberOfGenerators*SUNcAlgebra::NumberOfGenerators];\
SU_Nc_ADJOINT_FORMAT U2AdjzM[SUNcAlgebra::NumberOfGenerators*SUNcAlgebra::NumberOfGenerators];

#define COMPUTE_NEIGHBORING_ADJOINTS(x,y,z) \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,0),U0Adjx); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,1),U1Adjy); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,2),U2Adjz); \
\
SUNcGroup::Operations::GetAdjoint(U->Get(x-1,y,z,0),U0AdjxM); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y-1,z,1),U1AdjyM); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z-1,2),U2AdjzM);

////////////////////////////////////////////////////////
// BUFFERS FOR COMPUTATION OF AVERAGE FIELD STRENGTHS //
////////////////////////////////////////////////////////

// MAGNETIC //
#define SET_AVG_MAGNETIC_FIELD_STRENGTH_BUFFERS() \
\
SET_ELEMENTARY_PLAQUETTE_BUFFERS(); \
SET_NEIGHBORING_PLAQUETTE_BUFFERS(); \
SET_EXTENDED_PLAQUETTE_BUFFERS(); \
\
SU_Nc_FUNDAMENTAL_FORMAT UzxAvg[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UxyAvg[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UyzAvg[SUNcGroup::MatrixSize]; \
\
DOUBLE B0Loc[SUNcAlgebra::NumberOfGenerators];\
DOUBLE B1Loc[SUNcAlgebra::NumberOfGenerators];\
DOUBLE B2Loc[SUNcAlgebra::NumberOfGenerators];\
\
DOUBLE B0SqrLoc,B1SqrLoc,B2SqrLoc;

// ELECTRIC //
#define SET_AVG_ELECTRIC_FIELD_STRENGTH_BUFFERS() \
\
DOUBLE E0Loc[SUNcAlgebra::NumberOfGenerators];\
DOUBLE E1Loc[SUNcAlgebra::NumberOfGenerators];\
DOUBLE E2Loc[SUNcAlgebra::NumberOfGenerators];\
\
DOUBLE E0SqrLoc,E1SqrLoc,E2SqrLoc;

// ALL //
#define SET_AVG_FIELD_STRENGTH_BUFFERS() \
\
SET_AVG_MAGNETIC_FIELD_STRENGTH_BUFFERS();\
SET_AVG_ELECTRIC_FIELD_STRENGTH_BUFFERS();


/////////////////////////////////////////////////////////////
// COMPUTATION OF AVERAGE FIELD STRENGTHS AT POINT (X,Y,Z) //
/////////////////////////////////////////////////////////////

// MAGNETIC //
#define COMPUTE_AVG_MAGNETIC_FIELD_STRENGTH(x,y,z,U) \
\
COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z,U); \
COMPUTE_NEIGHBORING_PLAQUETTES(x,y,z,U); \
COMPUTE_EXTENDED_PLAQUETTES(x,y,z,U); \
\
SUNcGroup::AvgMatrix(Uzx,UMzx,UzMx,UMzMx,UzxAvg); \
SUNcGroup::AvgMatrix(Uxy,UMxy,UxMy,UMxMy,UxyAvg); \
SUNcGroup::AvgMatrix(Uyz,UMyz,UyMz,UMyMz,UyzAvg); \
\
SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UyzAvg,B0Loc); \
SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UzxAvg,B1Loc); \
SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UxyAvg,B2Loc); \
\
B0SqrLoc=DOUBLE(4.0)*SUNcGroup::Operations::ReTrIDMinusU(UyzAvg); \
B1SqrLoc=DOUBLE(4.0)*SUNcGroup::Operations::ReTrIDMinusU(UzxAvg); \
B2SqrLoc=DOUBLE(4.0)*SUNcGroup::Operations::ReTrIDMinusU(UxyAvg);


// ELECTRIC //
#define COMPUTE_AVG_ELECTRIC_FIELD_STRENGTH(x,y,z)\
\
SU_Nc_ALGEBRA_FORMAT E0Down[SUNcAlgebra::NumberOfGenerators]; \
SU_Nc_ALGEBRA_FORMAT E1Down[SUNcAlgebra::NumberOfGenerators]; \
SU_Nc_ALGEBRA_FORMAT E2Down[SUNcAlgebra::NumberOfGenerators]; \
\
SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x-1,y,z,0),E->Get(x-1,y,z,0,0),E0Down);\
SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x,y-1,z,1),E->Get(x,y-1,z,1,0),E1Down);\
SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x,y,z-1,2),E->Get(x,y,z-1,2,0),E2Down);\
\
for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){\
E0Loc[a]=DOUBLE(0.5)*(E->Get(x,y,z,0,a)[0]+E0Down[a]); \
E1Loc[a]=DOUBLE(0.5)*(E->Get(x,y,z,1,a)[0]+E1Down[a]); \
E2Loc[a]=DOUBLE(0.5)*(E->Get(x,y,z,2,a)[0]+E2Down[a]); \
} \
\
E0SqrLoc=0.0;   E1SqrLoc=0.0;   E2SqrLoc=0.0; \
\
for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){\
\
E0SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,0,a)[0])+SQR(E->Get(x-1,y,z,0,a)[0])); \
E1SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,1,a)[0])+SQR(E->Get(x,y-1,z,1,a)[0])); \
E2SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,2,a)[0])+SQR(E->Get(x,y,z-1,2,a)[0]));} \
\


/*// ELECTRIC //
 #define COMPUTE_AVG_ELECTRIC_FIELD_STRENGTH(x,y,z) \
 \
 for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){\
 E0Loc[a]=0.0;   E1Loc[a]=0.0;   E2Loc[a]=0.0;\
 }\
 \
 COMPUTE_NEIGHBORING_ADJOINTS(x,y,z); \
 \
 for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){\
 for(int b=0;b<SUNcAlgebra::NumberOfGenerators;b++){ \
 E0Loc[a]+=DOUBLE(0.5)*(DELTA(a,b)*E->Get(x,y,z,0,b)[0]+U0AdjxM[SUNcGroup::AdjIndex(b,a)]*E->Get(x-1,y,z,0,b)[0]); \
 E1Loc[a]+=DOUBLE(0.5)*(DELTA(a,b)*E->Get(x,y,z,1,b)[0]+U1AdjyM[SUNcGroup::AdjIndex(b,a)]*E->Get(x,y-1,z,1,b)[0]); \
 E2Loc[a]+=DOUBLE(0.5)*(DELTA(a,b)*E->Get(x,y,z,2,b)[0]+U2AdjzM[SUNcGroup::AdjIndex(b,a)]*E->Get(x,y,z-1,2,b)[0]); \
 } \
 }\
 \
 E0SqrLoc=0.0;   E1SqrLoc=0.0;   E2SqrLoc=0.0; \
 \
 for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){\
 \
 E0SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,0,a)[0])+SQR(E->Get(x-1,y,z,0,a)[0])); \
 E1SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,1,a)[0])+SQR(E->Get(x,y-1,z,1,a)[0])); \
 E2SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,2,a)[0])+SQR(E->Get(x,y,z-1,2,a)[0])); \
 \
 }*/

// ALL //
#define COMPUTE_AVG_FIELD_STRENGTH(x,y,z,U) \
\
COMPUTE_AVG_MAGNETIC_FIELD_STRENGTH(x,y,z,U);\
COMPUTE_AVG_ELECTRIC_FIELD_STRENGTH(x,y,z);
#endif

