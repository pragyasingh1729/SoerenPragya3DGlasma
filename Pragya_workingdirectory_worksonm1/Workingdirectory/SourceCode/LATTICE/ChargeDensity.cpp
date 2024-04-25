#ifndef __CHARGE_DENSITY_OBJECT__CPP__
#define __CHARGE_DENSITY_OBJECT__CPP__

#include <cmath>

///////////////////////////////////////////////////////////////////////////////
//                    WE DEFINE THE DIMENSIONLESS CURRENT VARIABLES AS       //
//                    J^{0}(x)=\psi^\dagger * t^a *\psi                     //
///////////////////////////////////////////////////////////////////////////////

class ChargeDensity{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
    
    DOUBLE a[Dimension];
    
    DOUBLE aCube;
    
    INT N[Dimension];

    SU_Nc_ALGEBRA_FORMAT *J;
    
    INT Index3D(INT x,INT y,INT z){
        
        return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*(z+1));
	   
    }

    INT Index(INT x,INT y,INT z,INT a){
        
        return a+SUNcAlgebra::NumberOfGenerators*Index3D(x,y,z);
        
    }
    
    SU_Nc_ALGEBRA_FORMAT* Get(INT x,INT y,INT z,INT a){
        
        return &J[Index(x,y,z,a)];
    }

    
    SU_Nc_ALGEBRA_FORMAT GetValue(INT x,INT y,INT z,INT a){
        
        return J[Index(x,y,z,a)];
    }
    
    void SynchronizeGhostCells(){
			   
    #if(BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)
	   //std::cerr << "#SYNCHRONIZING BOUNDARIES -- OPEN BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1,0),           this->Get(0,0,0,0)		      ,SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2],0),this->Get(0,0,Lattice::N[2]-1,0),SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
    #endif
			   
    #if(BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   std::cerr << "#SYNCHRNOZING BOUNDARIES -- PERIODIC BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1,0),           this->Get(0,0,Lattice::N[2]-1,0),SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2],0),this->Get(0,0,0,0)              ,SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
    #endif
			   
    }
    
    void SetZero(){
			   
	   #pragma omp parallel for schedule(static) collapse(3)
		for(INT z=0;z<=N[2]-1;z++){
		   for(INT y=0;y<=N[1]-1;y++){
			 for(INT x=0;x<=N[0]-1;x++){
				
			    for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				  
				  this->Get(x,y,z,a)[0]=DOUBLE(0.0);
				  
				}
			 }
		  }
	   }
	   this->SynchronizeGhostCells();
    }
			   

   ChargeDensity(INT Nx,INT Ny,INT Nz,DOUBLE ax,DOUBLE ay,DOUBLE az){
    
		this->Volume=Nx*Ny*Nz;
		
		this->N[0]=Nx;
		this->N[1]=Ny;
		this->N[2]=Nz;
		
		this->a[0]=ax;
		this->a[1]=ay;
		this->a[2]=az;
		
		this->aCube=a[0]*a[1]*a[2];
		
		// GHOST CELL SCHEME //
		this->J=(SU_Nc_ALGEBRA_FORMAT *) _mm_malloc(SUNcAlgebra::NumberOfGenerators*Nx*Ny*(Nz+2)*sizeof(SU_Nc_ALGEBRA_FORMAT),64);
			 
		this->SetZero();
		
      }
		
      ~ChargeDensity(){
		
       _mm_free(J);
    }
};
#endif
