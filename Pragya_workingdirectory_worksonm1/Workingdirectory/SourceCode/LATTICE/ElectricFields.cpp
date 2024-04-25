#ifndef __ELECTRICFIELDS__CPP__
#define __ELECTRICFIELDS__CPP__

////////////////////////////////////////////////////////////////////////////////////////////////
//                    WE DEFINE THE DIMENSIONLESS ELECTRIC FIELD VARIABLES AS                 //
//                 E^{\mu}(x)=g \sqrt{-g(x)} a^3/a_{mu}  g^{\mu\nu} \partial_{\tau} A_{\nu}   //
////////////////////////////////////////////////////////////////////////////////////////////////

class ElectricFields{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
    
    DOUBLE a[Dimension];
    
    DOUBLE aCube;
    
    INT N[Dimension];
        
    SU_Nc_ALGEBRA_FORMAT *E;
    
    INT Index3D(INT x,INT y,INT z){
        
        return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*(z+1));
        
    }
    
    INT Index(INT x,INT y,INT z,INT mu,INT a){
        
        return a+SUNcAlgebra::NumberOfGenerators*(mu+Dimension*Index3D(x,y,z));
        
    }
    
    SU_Nc_ALGEBRA_FORMAT* Get(INT x,INT y,INT z,INT mu,INT a){
        
        return &E[Index(x,y,z,mu,a)];
    }
    
    SU_Nc_ALGEBRA_FORMAT GetValue(INT x,INT y,INT z,INT mu,INT a){
        
        return E[Index(x,y,z,mu,a)];
    }

    void SynchronizeGhostCells(){
	   
        #if(BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)
	   //std::cerr << "#SYNCHRNIZING BOUNDARIES -- OPEN BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1,0,0),           this->Get(0,0,0,0,0)		       ,Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2],0,0),this->Get(0,0,Lattice::N[2]-1,0,0),Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
        #endif
	   
	   #if(BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   //std::cerr << "#SYNCHRNIZING BOUNDARIES -- PERIODIC BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1,0,0),           this->Get(0,0,Lattice::N[2]-1,0,0),Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2],0,0),this->Get(0,0,0,0,0)              ,Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   #endif
	   
    }
    
    void SynchronizeGhostCells(char D){
	   
	   if(D=='L'){
		  std::memcpy(this->Get(0,0,-1,0,0),           this->Get(0,0,0,0,0)		        ,Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   }
	   
	   if(D=='R'){
		  std::memcpy(this->Get(0,0,Lattice::N[2],0,0),this->Get(0,0,Lattice::N[2]-1,0,0)  ,Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_ALGEBRA_FORMAT));
	   }
	   
	   else
		 this->SynchronizeGhostCells();
		  
    }
    
    
    void SetZero(){
        
        #pragma omp parallel for schedule(static) collapse(3)
        for(INT z=0;z<=N[2]-1;z++){
            for(INT y=0;y<=N[1]-1;y++){
                for(INT x=0;x<=N[0]-1;x++){
                    
                    for(int mu=0;mu<Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            
                            this->Get(x,y,z,mu,a)[0]=DOUBLE(0.0);
                        
                        }
                    }
                    
                }
            }
        }
	   this->SynchronizeGhostCells();
    }
    
    
    ElectricFields(INT Nx,INT Ny,INT Nz,DOUBLE ax,DOUBLE ay,DOUBLE az){
        
        this->Volume=Nx*Ny*Nz;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        this->N[2]=Nz;
        
        this->a[0]=ax;
        this->a[1]=ay;
        this->a[2]=az;
        
        this->aCube=a[0]*a[1]*a[2];
        
        // GHOST CELL SCHEME //
        this->E=(SU_Nc_ALGEBRA_FORMAT *) _mm_malloc(Dimension*SUNcAlgebra::NumberOfGenerators*Nx*Ny*(Nz+2)*sizeof(SU_Nc_ALGEBRA_FORMAT),64);
        
        //this->E=(SU_Nc_ALGEBRA_FORMAT *) _mm_malloc(Dimension*SUNcAlgebra::NumberOfGenerators*Volume*sizeof(SU_Nc_ALGEBRA_FORMAT),64);
        
        this->SetZero();
        
    }
    
    ~ElectricFields(){
        
        _mm_free(E);
    }
};
/*
void Copy(ElectricFields *ECopy,ElectricFields *EOriginal){
    
    for(INT i=0;i<EOriginal->Dimension;i++){
        
        ECopy->a[i]=EOriginal->a[i];
        
        ECopy->N[i]=EOriginal->N[i];
        
    }
    
    ECopy->Volume=EOriginal->Volume;
    
    ECopy->aCube=EOriginal->aCube;
    
    std::memcpy(ECopy->Get(0,0,0,0,0),EOriginal->Get(0,0,0,0,0),Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*Lattice::Volume*sizeof(SU_Nc_ALGEBRA_FORMAT));
    
    
}
*/
/////////////////////////
//BLOCK ELECTRIC FIELDS//
/////////////////////////

void Block(ElectricFields **EOld){
    
    // CREATE NEW OBJECT //
    ElectricFields *ENew=new ElectricFields((*EOld)->N[0]/2,(*EOld)->N[1]/2,(*EOld)->N[2]/2,2*(*EOld)->a[0],2*(*EOld)->a[1],2*(*EOld)->a[2]);
    
    // SET ELECTRIC FIELDS TO ZERO //
    ENew->SetZero();
    
    // DELETE OLD OBJECT //
    delete *EOld;
    
    // SET POINTER TO BLOCKED FIELDS //
    *EOld=ENew;
    
}

#endif

