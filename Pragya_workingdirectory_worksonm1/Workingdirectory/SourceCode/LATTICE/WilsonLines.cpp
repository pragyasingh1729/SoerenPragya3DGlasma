#ifndef __WILSONLINES__CPP__
#define __WILSONLINES__CPP__

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WE DEFINE WILSON LINES AS                                               //
//                          V^\DAGGER(x_plus,x)= P exp[ig 1/RHO(x_plus,x)]                              //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class WilsonLines{
    
public:

    INT Volume;
    
    const static INT Dimension=3;
            
    DOUBLE a[Dimension];
    
    DOUBLE aCube;
    
    INT N[Dimension];
        
    SU_Nc_FUNDAMENTAL_FORMAT *V;

    INT Index3D(INT x,INT y,INT z){

         return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*(z+1));
    }

    INT Index(INT x,INT y, INT z){

        return SUNcGroup::MatrixSize*(Index3D(x,y,z));
    }

    SU_Nc_FUNDAMENTAL_FORMAT* Get(INT x,INT y,INT z){
        
        return &(V[Index(x,y,z)]);
    }
    
    void SynchronizeGhostCells(){

    #if(BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)
	   //std::cerr << "#SYNCHRNIZING BOUNDARIES -- OPEN BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1),            this->Get(0,0,0)			      ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2]), this->Get(0,0,Lattice::N[2]-1)   ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
    #endif

    #if(BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   std::cerr << "#SYNCHRNIZING BOUNDARIES -- PERIODIC BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1),            this->Get(0,0,Lattice::N[2]-1)   ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2]), this->Get(0,0,0)                 ,SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
    #endif

    }
  
    void SetIdentity(){
        
        #pragma omp parallel for schedule(static) collapse(3)
        for(INT z=0;z<=N[2]-1;z++){
            for(INT y=0;y<=N[1]-1;y++){
                for(INT x=0;x<=N[0]-1;x++){
                    for(int mu=0;mu<Dimension;mu++){

                        COPY_SUNcMatrix(this->Get(x,y,z),SUNcGroup::UnitMatrix);
                    
                    }
                }
            }
        }
	  
	  this->SynchronizeGhostCells();
    }

    WilsonLines(INT Nx,INT Ny,INT Nz,DOUBLE ax,DOUBLE ay,DOUBLE az){
        
        this->Volume=Nx*Ny*Nz;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        this->N[2]=Nz;
        
        this->a[0]=ax;
        this->a[1]=ay;
        this->a[2]=az;
        
        this->aCube=a[0]*a[1]*a[2];
        
	   this->V=(SU_Nc_FUNDAMENTAL_FORMAT *) _mm_malloc(Lattice::Dimension*SUNcGroup::MatrixSize*Nx*Ny*(Nz+2)*sizeof(SU_Nc_FUNDAMENTAL_FORMAT),64);
        
        this->SetIdentity();

    }

    ~WilsonLines(){
        
        _mm_free(V);
    }
    
    
};

#endif



// Part coming from Charge density to be added later line 92
/*
 void Print(){
 
 for(INT z=0;z<=N[2]-1;z++){
 for(INT y=0;y<=N[1]-1; y++){
 for (INT x=0; x<=N[0]-1;x++){
 
 std::cout << x << " " << y << " " << z;
 
 
 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
 
 std::cout << " "<< this->GetValue(x,y,z,a);
 
 }
 
 
 std::cout << std::endl;
 }
 
 std::cout << std::endl;
 
 }
 
 std::cout << std::endl;
 
 }
 
 
 }
 
 void Output(std::string fname){
 
 std::ofstream OutStream;
 OutStream.open(fname.c_str());
 
 for(INT z=0;z<=N[2]-1; z++){
 for(INT y=0;y<=N[1]-1; y++){
 for (INT x=0; x<=N[0]-1;x++){
 
 OutStream << x << " " << y << " " << z;
 
 
 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
 
 OutStream << " "<< this->GetValue(x,y,z,a);
 
 }
 
 
 OutStream<<std::endl;
 }
 
 OutStream<<std::endl;
 
 }
 
 OutStream<<std::endl;
 
 }
 
 OutStream.close();
 }
 
 */
