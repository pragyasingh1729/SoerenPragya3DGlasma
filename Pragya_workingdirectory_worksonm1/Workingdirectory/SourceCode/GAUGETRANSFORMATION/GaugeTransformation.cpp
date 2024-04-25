#ifndef __GAUGETRANSFORMATION__CPP__
#define __GAUGETRANSFORMATION__CPP__

////////////////////////////////////////////////
//GAUGE TRANSFORMED LINKS AND ELECTRIC FIELDS //
////////////////////////////////////////////////

namespace GaugeFixedVariables{
    
    GaugeLinks *U;
    ElectricFields *E;
    
    //INITIALIZE VARIABLES
    void Init(){
               
        U=new GaugeLinks(SUNcGaugeLinks::U->N[0],SUNcGaugeLinks::U->N[1],SUNcGaugeLinks::U->N[2],SUNcGaugeLinks::U->a[0],SUNcGaugeLinks::U->a[1],SUNcGaugeLinks::U->a[2]);
        E=new ElectricFields(SUNcElectricFields::E->N[0],SUNcElectricFields::E->N[1],SUNcElectricFields::E->N[2],SUNcElectricFields::E->a[0],SUNcElectricFields::E->a[1],SUNcElectricFields::E->a[2]);

    }
    
}

/////////////////////////////////////////////////////////////////////////
//GAUGE TRANSFORMATION AND OPERATIONS TO PERFORM GAUGE TRANSFORMATIONS //
/////////////////////////////////////////////////////////////////////////

namespace GaugeTransformation{
    
    GaugeTransformations *G;

    ///////////////////
    //INITIALIZATION //
    ///////////////////
    
    void Init(){
        
        G=new GaugeTransformations(SUNcGaugeLinks::U->N[0],SUNcGaugeLinks::U->N[1],SUNcGaugeLinks::U->N[2],SUNcGaugeLinks::U->a[0],SUNcGaugeLinks::U->a[1],SUNcGaugeLinks::U->a[2]);

    }
    
    void SetIdentity(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        #pragma omp parallel for
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    COPY_SUNcMatrix(G->Get(x,y,z),SUNcGroup::UnitMatrix);
                    
                }
            }
        }// END PARALLEL FOR
	   G->SynchronizeGhostCells();
    }
    
    void SetIdentity(){
        
        SetIdentity(0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1);
        
    }
    
    void SetRandom(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    RandomNumberGenerator::SUNcMatrix(DOUBLE(2.0)*M_SQRT2*M_PI,G->Get(x,y,z));
                    
                }
            }
        }
         G->SynchronizeGhostCells();
    }
    
    void SetRandom(){
        
        SetRandom(0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1);
    
    }
    
    void Save(std::string fname){
        Save(fname,G);
    }
    
    //////////////////////////
    //GAUGE TRANSFORMATIONS //
    //////////////////////////
    
    namespace Operations{
        
        //CREATE A COPY OF THE DYNAMICAL FIELDS
        void CopyFields(){
            
            std::memcpy(GaugeFixedVariables::U->Get(0,0,0,0),SUNcGaugeLinks::U->Get(0,0,0,0),Lattice::Dimension*SUNcGroup::MatrixSize*SUNcGaugeLinks::U->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
            std::memcpy(GaugeFixedVariables::E->Get(0,0,0,0,0),SUNcElectricFields::E->Get(0,0,0,0,0),Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*SUNcElectricFields::E->Volume*sizeof(SU_Nc_ALGEBRA_FORMAT));
        }
        
        //CREATE A COPY OF THE DYNAMICAL FIELDS
        void SaveFields(){
            
            std::memcpy(SUNcGaugeLinks::U->Get(0,0,0,0),GaugeFixedVariables::U->Get(0,0,0,0),Lattice::Dimension*SUNcGroup::MatrixSize*SUNcGaugeLinks::U->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
            std::memcpy(SUNcElectricFields::E->Get(0,0,0,0,0),GaugeFixedVariables::E->Get(0,0,0,0,0),Lattice::Dimension*SUNcAlgebra::NumberOfGenerators*SUNcElectricFields::E->Volume*sizeof(SU_Nc_ALGEBRA_FORMAT));
            
		  G->SynchronizeGhostCells();

        }
        

        //PERFORM GAUGE TRANSFORMATION OF GAUGE LINKS
        void GaugeTransformLinks(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
            
            #pragma omp parallel for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        PERFORM_GAUGE_LINK_TRANSFORMATION(x,y,z)
		
                    }
                }
            }// END PARALLEL FOR
		   G->SynchronizeGhostCells();
        }
        
        void GaugeTransformLinks(GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
            
            GaugeTransformLinks(0,UOld->N[0]-1,0,UOld->N[1]-1,0,UOld->N[2]-1,UOld,UNew,G);
            
        }
        
        void GaugeTransformLinks(){
            
            GaugeTransformLinks(SUNcGaugeLinks::U,GaugeFixedVariables::U,GaugeTransformation::G);
            
        }
        
        //PERFORM GAUGE TRANSFORMATION OF ELECTRIC FIELDS
        void GaugeTransformElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ElectricFields *EOld,ElectricFields *ENew,GaugeTransformations *G){
            
                #pragma omp parallel for
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            PERFORM_ELECTRIC_FIELD_TRANSFORMATION(x,y,z)
                            
                        }
                    }
                }
             G->SynchronizeGhostCells();
        }
        
        void GaugeTransformElectricFields(ElectricFields *EOld,ElectricFields *ENew,GaugeTransformations *G){
            
            GaugeTransformElectricFields(0,EOld->N[0]-1,0,EOld->N[1]-1,0,EOld->N[2]-1,EOld,ENew,G);
            
        }
        
        void GaugeTransformElectricFields(){
            
            GaugeTransformElectricFields(SUNcElectricFields::E,GaugeFixedVariables::E,GaugeTransformation::G);
            
        }
        
        
    }
	
}

#endif
