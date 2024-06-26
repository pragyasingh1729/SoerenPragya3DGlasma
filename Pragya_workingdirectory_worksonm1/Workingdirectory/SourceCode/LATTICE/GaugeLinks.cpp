#ifndef __GAUGELINKS__CPP__
#define __GAUGELINKS__CPP__

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WE DEFINE THE LATTICE GAUGE LINKS AS                                    //
//                          U^{\dagger}_{mu}(x)= P exp[ -ig a_{\mu} A_{\mu}(x)]                         //
//                               WITH DERIVATIVES ACTING AS                                             //
//   \delta U_{mu}(x)/ \delta A_{nu}^{a}(y)= -ig a_{\mu}t^{a} U_{mu}(x) \delta_{x,y} \delta_{\mu\nu}    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class GaugeLinks{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
    
    DOUBLE a[Dimension];
    
    DOUBLE aCube;
    
    INT N[Dimension];
    
    SU_Nc_FUNDAMENTAL_FORMAT *U;
    
    INT Index3D(INT x,INT y,INT z){
	   
	   if(z<-1 || z>this->N[2]+1){
		  std::cerr << "#CRITICAL ERRROR -- z " << z << "OUT OF RANGE " << std::endl;
		  exit(0);
	   }
	   
	   //return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*MOD(z,N[2]));
	   return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*(z+1));
    }
    
    INT Index(INT x,INT y,INT z,INT mu){
	   
	   return SUNcGroup::MatrixSize*(mu+Dimension*Index3D(x,y,z));
    }
    
    // void GetPosition(INT Index,INT &x,INT &y,INT &z){
    
    //     z=Index/(N[0]*N[1]); y=(Index-z*(N[0]*N[1]))/(N[0]); x=(Index-N[0]*(y+N[1]*z));
    
    // }
    
    SU_Nc_FUNDAMENTAL_FORMAT* Get(INT x,INT y,INT z,INT mu){
	   
	   return &(U[Index(x,y,z,mu)]);
    }
    
    SU_Nc_FUNDAMENTAL_FORMAT* Get(INT xyzPos,INT mu){
	   
	   return &(U[SUNcGroup::MatrixSize*(mu+Dimension*(xyzPos))]);
    }
    
    void SynchronizeGhostCells(){
	   
	#if(BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)

	   std::memcpy(this->Get(0,0,-1,0),           this->Get(0,0,0,0)			  ,Lattice::Dimension*SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2],0),this->Get(0,0,Lattice::N[2]-1,0),Lattice::Dimension*SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	#endif
	   
	#if(BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   std::cerr << "#SYNCHRNIZING BOUNDARIES -- PERIODIC BOUNDARY CONDITIONS" << std::endl;
	   std::memcpy(this->Get(0,0,-1,0),           this->Get(0,0,Lattice::N[2]-1,0),Lattice::Dimension*SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   std::memcpy(this->Get(0,0,Lattice::N[2],0),this->Get(0,0,0,0)              ,Lattice::Dimension*SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	#endif
	   
    }
    
    void SynchronizeGhostCells(char D){

	   if(D=='L'){
		  std::memcpy(this->Get(0,0,-1,0),           this->Get(0,0,0,0)			   ,Lattice::Dimension*SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   }
	   
	   if(D=='R'){
		  std::memcpy(this->Get(0,0,Lattice::N[2],0),this->Get(0,0,Lattice::N[2]-1,0),Lattice::Dimension*SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
	   }
	   
	  else
		 this->SynchronizeGhostCells();
    }
    
    
    void SetIdentity(){
	   
    #pragma omp parallel for schedule(static) collapse(3)
	   for(INT z=0;z<=N[2]-1;z++){
		  for(INT y=0;y<=N[1]-1;y++){
			 for(INT x=0;x<=N[0]-1;x++){
				for(int mu=0;mu<Dimension;mu++){
				    
				    COPY_SUNcMatrix(this->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
				    
				}
			 }
		  }
	   }
	   
	   this->SynchronizeGhostCells();
    }
    
    
    
    GaugeLinks(INT Nx,INT Ny,INT Nz,DOUBLE ax,DOUBLE ay,DOUBLE az){
	   
	   this->Volume=Nx*Ny*Nz;
	   
	   this->N[0]=Nx;
	   this->N[1]=Ny;
	   this->N[2]=Nz;
	   
	   this->a[0]=ax;
	   this->a[1]=ay;
	   this->a[2]=az;
	   
	   this->aCube=a[0]*a[1]*a[2];
	   
	   this->U=(SU_Nc_FUNDAMENTAL_FORMAT *) _mm_malloc(Dimension*SUNcGroup::MatrixSize*Nx*Ny*(Nz+2)*sizeof(SU_Nc_FUNDAMENTAL_FORMAT),64);
	   this->SetIdentity();
	   
    }
    
    ~GaugeLinks(){
	   
	   _mm_free(U);
    }
    
    
};


//void Copy(GaugeLinks *UCopy,GaugeLinks *UOriginal){
//
//    for(INT i=0;i<UOriginal->Dimension;i++){
//
//        UCopy->a[i]=UOriginal->a[i];
//
//        UCopy->N[i]=UOriginal->N[i];
//
//    }
//    //UCopy->(Lattice::N[0]*Lattice::N[1]*(Lattice::N[2]+2))=UOriginal->(Lattice::N[0]*Lattice::N[1]*(Lattice::N[2]+2));
//    //UCopy->Volume=UOriginal->Volume;
//    UCopy->TotalVolume=UOriginal->TotalVolume;
//    UCopy->aCube=UOriginal->aCube;
//
//    // TODO: NEED TO COPY GHOST CELLS ALSO REPLACE VOLUME BY ARRAY SIZE //
//    //std::memcpy(UCopy->Get(0,0,0,0),UOriginal->Get(0,0,0,0),UOriginal->Dimension*SUNcGroup::MatrixSize*UOriginal->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
//    std::memcpy(UCopy->Get(0,0,0,0),UOriginal->Get(0,0,0,0),UOriginal->Dimension*SUNcGroup::MatrixSize*UOriginal->TotalVolume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
//
//}


/*
 void CopyWithBlock(GaugeLinks **UOld,GaugeLinks *UNew){
 
 // COMPUTE BLOCKED GAUGE LINKS //
 for(INT z=0;z<=UNew->N[2]-1;z++){
 for(INT y=0;y<=UNew->N[1]-1;y++){
 for(INT x=0;x<=UNew->N[0]-1;x++){
 
 SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,0),UNew->Get(x,y,z,0));
 SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,1),UNew->Get(x,y,z,1));
 SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,2),UNew->Get(x,y,z,2));
 
 }
 }
 }
 
 }
 
 void Block(GaugeLinks **UOld){
 
 // CREATE NEW OBJECTS //
 GaugeLinks *UNew=new GaugeLinks((*UOld)->N[0]/2,(*UOld)->N[1]/2,(*UOld)->N[2]/2,2*(*UOld)->a[0],2*(*UOld)->a[1],2*(*UOld)->a[2]);
 
 // COMPUTE BLOCKED GAUGE LINKS //
 for(INT z=0;z<=UNew->N[2]-1;z++){
 for(INT y=0;y<=UNew->N[1]-1;y++){
 for(INT x=0;x<=UNew->N[0]-1;x++){
 
 SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,0),UNew->Get(x,y,z,0));
 SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,1),UNew->Get(x,y,z,1));
 SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,2),UNew->Get(x,y,z,2));
 
 }
 }
 }
 
 // DELETE OLD OBJECTS //
 delete *UOld;
 
 // SET POINTER TO BLOCKED FIELDS //
 *UOld=UNew;
 
 }*/

#endif

