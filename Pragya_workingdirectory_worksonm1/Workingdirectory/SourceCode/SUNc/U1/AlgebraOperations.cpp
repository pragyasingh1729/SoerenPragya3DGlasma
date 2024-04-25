#ifndef __SU_Nc_ALGEBRA_OPERATIONS__
#define __SU_Nc_ALGEBRA_OPERATIONS__

namespace SUNcAlgebra{
    
    //U(1) GENERATORS ARE t=1/sqrt(2)
    static const int NumberOfGenerators=1;
    
    static const COMPLEX GeneratorMatrix[Nc*Nc][NumberOfGenerators]={{1.0/M_SQRT2}};
    
    //COMPUTES MATRIX REPRESENTATION  it^{a}alpha^a ///
    void GetMatrixFormIAlpha(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *alphaMatrix){
        
        alphaMatrix[0]=c*COMPLEX(0.0,0.5)*DOUBLE(M_SQRT2)*alpha[0];
        
    }
    
    //BASIC OPERATIONS
    namespace Operations{
        
        //COMPUTES exp( c i t^{a} alpha^{a}) BY DIAGONALIZING THE MATRIX
        void MatrixIExp(DOUBLE c,DOUBLE *alpha,SU_Nc_FUNDAMENTAL_FORMAT *ExpIAlpha){
            
            ExpIAlpha[0]=COMPLEX (cos(DOUBLE(0.5)*c*M_SQRT2*alpha[0]),sin(DOUBLE(0.5)*c*M_SQRT2*alpha[0]));
            
        }
        
        
        void MatrixILog(DOUBLE c,SU_Nc_FUNDAMENTAL_FORMAT *U, SU_Nc_ALGEBRA_FORMAT *LogU){
            
            LogU[0]=M_SQRT2/c*arg(U[0]);
            
        }
        
        // COMPUTES t^{a}.U  //
        void GeneratorProduct(COMPLEX c,SU_Nc_FUNDAMENTAL_FORMAT *U, COMPLEX (&Utau)[Nc*Nc][NumberOfGenerators]){
            Utau[0][0]=c*U[0]/M_SQRT2;
            
        }
        
        void AdjointMultiplication(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ALGEBRA_FORMAT *E,SU_Nc_ALGEBRA_FORMAT *ENew){
            
            ENew[0]=E[0];
        }
        
        void InverseAdjointMultiplication(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ALGEBRA_FORMAT *E,SU_Nc_ALGEBRA_FORMAT *ENew){
            
            ENew[0]=E[0]; 
        }
        
        
        
        
        
    }
    
}

#endif
