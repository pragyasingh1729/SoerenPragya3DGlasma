
#ifndef __GAUGETRANSFORMATIONS__CPP__
#define __GAUGETRANSFORMATIONS__CPP__

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WE DEFINE THE LATTICE GAUGE LINKS AS                                    //
//                          U^{\dagger}_{mu}(x)= P exp[ -ig a_{\mu} A_{\mu}(x)]                         //
//                               WITH DERIVATIVES ACTING AS                                             //
//   \delta U_{mu}(x)/ \delta A_{nu}^{a}(y)= -ig a_{\mu}t^{a} U_{mu}(x) \delta_{x,y} \delta_{\mu\nu}    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class GaugeTransformations{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
    
    DOUBLE a[Dimension];
    
    INT N[Dimension];
        
    SU_Nc_FUNDAMENTAL_FORMAT *G;
    
    INT Index3D(INT x,INT y,INT z){
        
        //return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*MOD(z,N[2]));
	   return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*(z+1));

    }
    
    INT Index(INT x,INT y,INT z){

        return SUNcGroup::MatrixSize*(Index3D(x,y,z));
    
    }
    
    SU_Nc_FUNDAMENTAL_FORMAT* Get(INT x,INT y,INT z){
        
        return &(G[Index(x,y,z)]);
    }
    
    void SynchronizeGhostCells(){
	   
#if(BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)
	   
	   std::memcpy(this->Get(0,0,-1),           this->Get(0,0,0)			      ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2]),this->Get(0,0,Lattice::N[2]-1)    ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
#endif
	   
#if(BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   std::cerr << "#SYNCHRNIZING BOUNDARIES -- PERIODIC BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1),           this->Get(0,0,Lattice::N[2]-1)   ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2]),this->Get(0,0,0)                 ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
#endif
	   
    }
    
    void SetIdentity(){
        
        #pragma omp parallel for schedule(static) collapse(3)
        for(INT z=0;z<=N[2]-1;z++){
            for(INT y=0;y<=N[1]-1;y++){
                for(INT x=0;x<=N[0]-1;x++){
                    
                    COPY_SUNcMatrix(this->Get(x,y,z),SUNcGroup::UnitMatrix);
                    
                }
            }
        }
	   
	   this->SynchronizeGhostCells();
    }
    
    GaugeTransformations(INT Nx,INT Ny,INT Nz,DOUBLE ax,DOUBLE ay,DOUBLE az){
        
        this->Volume=Nx*Ny*Nz;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        this->N[2]=Nz;
        
        this->a[0]=ax;
        this->a[1]=ay;
        this->a[2]=az;
        
        this->G=(SU_Nc_FUNDAMENTAL_FORMAT *) _mm_malloc(SUNcGroup::MatrixSize*Nx*Ny*(Nz+2)*sizeof(SU_Nc_FUNDAMENTAL_FORMAT),64);
        
        this->SetIdentity();
        
    }
    
    ~GaugeTransformations(){
        
        _mm_free(G);
    }
    
    
};

void Save(std::string fname,GaugeTransformations *GaugeTrafo){
    
    std::ofstream OutStream;
    
    OutStream.open(fname.c_str());
    
    for(INT z=0;z<GaugeTrafo->N[2];z++){
        for(INT y=0;y<GaugeTrafo->N[1];y++){
            for(INT x=0;x<GaugeTrafo->N[0];x++){
                
                OutStream << x << " " << y << " " << z << " " << SUNcGroup::IO::MatrixToString(GaugeTrafo->Get(x,y,z)) << std::endl;
                
            }
        }
    }
    
    OutStream.close();    
}

#endif
