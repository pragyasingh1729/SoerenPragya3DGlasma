#ifndef __SU_NC_ALGEBRA_OPERATIONS__
#define __SU_NC_ALGEBRA_OPERATIONS__

namespace SUNcAlgebra{
    
    //SU(3) GENERATORS ARE t^{a}=lambda_{a}/2 WHERE lambda_{a} (a=1,..,8) ARE THE GELL-MANN MATRICES
    static const int NumberOfGenerators=8;
    
    //COMPUTES MATRIX REPRESENTATION  it^{a} alpha^a ///
    void GetMatrixFormIAlpha(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *alphaMatrix){
        
        alphaMatrix[SUNcGroup::MatrixIndex(0,0)]= c*COMPLEX(0.0,0.5)*(alpha[2]+alpha[7]/M_SQRT3);
        alphaMatrix[SUNcGroup::MatrixIndex(1,1)]=-c*COMPLEX(0.0,0.5)*(alpha[2]-alpha[7]/M_SQRT3);
        alphaMatrix[SUNcGroup::MatrixIndex(2,2)]=-c*COMPLEX(0.0,alpha[7]/M_SQRT3);
        
        alphaMatrix[SUNcGroup::MatrixIndex(1,0)]=c*COMPLEX(0.0,0.5)*COMPLEX(alpha[0],alpha[1]);
        alphaMatrix[SUNcGroup::MatrixIndex(2,0)]=c*COMPLEX(0.0,0.5)*COMPLEX(alpha[3],alpha[4]);
        alphaMatrix[SUNcGroup::MatrixIndex(2,1)]=c*COMPLEX(0.0,0.5)*COMPLEX(alpha[5],alpha[6]);
        
        alphaMatrix[SUNcGroup::MatrixIndex(0,1)]=conj(alphaMatrix[SUNcGroup::MatrixIndex(1,0)]);
        alphaMatrix[SUNcGroup::MatrixIndex(0,2)]=conj(alphaMatrix[SUNcGroup::MatrixIndex(2,0)]);
        alphaMatrix[SUNcGroup::MatrixIndex(1,2)]=conj(alphaMatrix[SUNcGroup::MatrixIndex(2,1)]);
        
        
    }
        
    //BASIC OPERATIONS
    namespace Operations{    
        
        //COMPUTES MATRIX REPRESENTATION  t^{a} alpha^a ///
        void MatrixForm(SU_Nc_ALGEBRA_FORMAT *alpha,COMPLEX A[3][3]){
            
            A[0][0]= DOUBLE(0.5)*(alpha[2]+alpha[7]/M_SQRT3);
            A[1][1]=-DOUBLE(0.5)*(alpha[2]-alpha[7]/M_SQRT3);
            A[2][2]=-alpha[7]/M_SQRT3;
            
            A[1][0]=DOUBLE(0.5)*COMPLEX(alpha[0],alpha[1]);
            A[2][0]=DOUBLE(0.5)*COMPLEX(alpha[3],alpha[4]);
            A[2][1]=DOUBLE(0.5)*COMPLEX(alpha[5],alpha[6]);
            
            A[0][1]=conj(A[1][0]);
            A[0][2]=conj(A[2][0]);
            A[1][2]=conj(A[2][1]);
            
        }
        
        //COMPUTES exp( c i t^{a} alpha^{a}) BY DIAGONALIZING THE MATRIX
        void MatrixIExp(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *ExpIAlpha){
            
            //GET MATRIX REPRESENTATION
            COMPLEX A[3][3];
            MatrixForm(alpha,A);
            
            //COMPUTE EIGENSYSTEM
            COMPLEX U[3][3];
            DOUBLE eV[3];
            
            int SUCCESS=Diagonalization3x3::Eigensystem(A,U,eV);
            
            if(SUCCESS==-1){
                std::cerr << "DIAGONALIZAITON FAILED!" << std::endl;
                exit(0);
            }
            
            //COMPUTE MATRIX EXPONENTIAL
            
            //SET TO ZERO
            for(int i=0;i<3;i++){
                for(int k=0;k<3;k++){
                    ExpIAlpha[SUNcGroup::MatrixIndex(i,k)]=DOUBLE(0.0);
                }
            }
            
            //PERFORM EXPONENTIATION
            for(int j=0;j<3;j++){
                
                //COMPUTE EXPONENTIAL OF I TIMES THE EIGENVALUE TIMES THE CONSTANT FACTOR c
                COMPLEX ExpIEv=COMPLEX(cos(c*eV[j]),sin(c*eV[j]));
                
                //CONTRACT WITH MATRIX
                for(int i=0;i<3;i++){
                    for(int k=0;k<3;k++){
                        ExpIAlpha[SUNcGroup::MatrixIndex(i,k)]+=U[i][j]*ExpIEv*conj(U[k][j]);
                    }
                }
                
            }
            
        }
        
    }
    
}

#endif
