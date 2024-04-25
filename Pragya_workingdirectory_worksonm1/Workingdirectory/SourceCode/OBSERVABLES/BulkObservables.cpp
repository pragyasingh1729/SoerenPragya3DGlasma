#ifndef __BULKOBSERVABLES__CPP__
#define __BULKOBSERVABLES__CPP__

namespace Observables {
        
    namespace Bulk{
        
        // ELECTRIC AND MAGNETIC FIELD STRENGTH SQUARED
        DOUBLE ESqr[Lattice::Dimension];
        DOUBLE BSqr[Lattice::Dimension];
        
        // SOURCE TERM
        DOUBLE EDotJ[Lattice::Dimension];
        
        //CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
        DOUBLE cB[Lattice::Dimension];
        DOUBLE cE[Lattice::Dimension];
        DOUBLE cS[Lattice::Dimension];
        
        
        //SET CONSTANTS FOR ENERGY MOMENTUM TENSOR
        void SetConstants(GaugeLinks *U){

            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                cB[mu]=SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                cE[mu]=SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                cS[mu]=(U->a[mu]*(Lattice::aScale)/SQR(U->aCube));
            }
            
        }
        
        
        // DIAGONAL COMPONENTS OF ENERGY MOMENTUM TENSOR FOR ORIGINAL GaugeLinks
        DOUBLE T00(){
            return 0.5*(cE[0]*ESqr[0]+cE[1]*ESqr[1]+cE[2]*ESqr[2]+cB[0]*BSqr[0]+cB[1]*BSqr[1]+cB[2]*BSqr[2]);
        }
        
        DOUBLE TXX(){
            
            return 0.5*((cE[1]*ESqr[1]+cE[2]*ESqr[2]+cB[1]*BSqr[1]+cB[2]*BSqr[2])-(cE[0]*ESqr[0]+cB[0]*BSqr[0]));
        }
        
        DOUBLE TYY(){
            return 0.5*((cE[0]*ESqr[0]+cE[2]*ESqr[2]+cB[0]*BSqr[0]+cB[2]*BSqr[2])-(cE[1]*ESqr[1]+cB[1]*BSqr[1]));
        }
        
        DOUBLE TZZ(){
            return 0.5*((cE[0]*ESqr[0]+cE[1]*ESqr[1]+cB[0]*BSqr[0]+cB[1]*BSqr[1])-(cE[2]*ESqr[2]+cB[2]*BSqr[2]));
        }
        
        DOUBLE ELECTRIC(){
            return 0.5*(cE[0]*ESqr[0]+cE[1]*ESqr[1]+cE[2]*ESqr[2]);
        }
        
        DOUBLE MAGNETIC(){
            return 0.5*(cB[0]*BSqr[0]+cB[1]*BSqr[1]+cB[2]*BSqr[2]);
        }
        
        DOUBLE SOURCE(){
            return cS[0]*EDotJ[0]+cS[1]*EDotJ[1]+cS[2]*EDotJ[2];
        }
        
        
        void UpdateMagnetic(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U){
            
            //SET OBSERVABLES
            DOUBLE B0Sqr=0.0; DOUBLE B1Sqr=0.0; DOUBLE B2Sqr=0.0;
            
            #pragma omp parallel
            {
                //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
                SET_ELEMENTARY_PLAQUETTE_BUFFERS();
                
                //UPDATE ALL GAUGE LINKS
                #pragma omp for reduction( + : B0Sqr,B1Sqr,B2Sqr)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            //COMPUTE ELEMENTARY PLAQUETTES AT THIS POINT
                            COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z,U);
                            
                            //ADD SQUARED FIELD STRENGTH TO LOCAL OBSERVABLES
                            B0Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uyz);
                            B1Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uzx);
                            B2Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uxy);
                            
                        }
                    }
                }
                
            } // END PARALLEL
		
//		    BSqr[0]=B0Sqr/U->Volume; BSqr[1]=B1Sqr/U->Volume; BSqr[2]=B2Sqr/U->Volume;
            BSqr[0]=B0Sqr; BSqr[1]=B1Sqr; BSqr[2]=B2Sqr;
        }
        
        void UpdateMagnetic(GaugeLinks *U){
            
            //SET CONSTANTS TO COMPUTE ENERGY MOMENTUM TENSOR
            SetConstants(U);
            
            UpdateMagnetic(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);

        }
        
        void UpdateElectric(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ElectricFields *E){
            
            //SET OBSERVABLES
            DOUBLE E0Sqr=0.0; DOUBLE E1Sqr=0.0; DOUBLE E2Sqr=0.0;
            
            #pragma omp parallel
            {
                        
                //UPDATE ALL GAUGE LINKS
                #pragma omp for reduction( + : E0Sqr,E1Sqr,E2Sqr)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
					   //SUM ALL COLOR COMPONENTS AND ADD TO SQUARED FIELD STRENGTH
					   		for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						  
							  E0Sqr+=SQR(E->Get(x,y,z,0,a)[0]);
							  E1Sqr+=SQR(E->Get(x,y,z,1,a)[0]);
							  E2Sqr+=SQR(E->Get(x,y,z,2,a)[0]);
					    	}
				     	}
					}
                }
                
            } // END PARALLEL
		 
//		  ESqr[0]=E0Sqr/E->Volume; ESqr[1]=E1Sqr/E->Volume; ESqr[2]=E2Sqr/E->Volume;
            ESqr[0]=E0Sqr; ESqr[1]=E1Sqr; ESqr[2]=E2Sqr;
            
        }
        
        void UpdateSource(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ElectricFields *E,ChargeDensity *Jp,ChargeDensity *Jm){
            
            //SET OBSERVABLES
            DOUBLE E0DotJ0=0.0;
            DOUBLE E1DotJ1=0.0;
            DOUBLE E2DotJ2=0.0;

                //UPDATE ALL GAUGE LINKS
                #pragma omp parallel for reduction( + : E0DotJ0,E1DotJ1,E2DotJ2)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            //SUM ALL COLOR COMPONENTS AND ADD TO SQUARED FIELD STRENGTH
                            for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                                
                                E0DotJ0+=0.0;
                                E1DotJ1+=0.0;
                                E2DotJ2+=E->Get(x,y,z,2,a)[0]*(Jp->Get(x,y,z,a)[0]-Jm->Get(x,y,z,a)[0])/M_SQRT2;
                                
                            }
                            
                        }
                    }
                }
            
//            EDotJ[0]=E0DotJ0/E->Volume; EDotJ[1]=E1DotJ1/E->Volume;  EDotJ[2]=E2DotJ2/E->Volume;
		    EDotJ[0]=E0DotJ0; EDotJ[1]=E1DotJ1;  EDotJ[2]=E2DotJ2;

        }
        
        
        void Update(GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp,ChargeDensity *Jm){
            
            //SET CONSTANTS TO COMPUTE ENERGY MOMENTUM TENSOR
            SetConstants(U);
            
            //UPDATE INDIVIDUAL CONTRIBUTIONS
            UpdateMagnetic(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);
            UpdateElectric(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,E);
            UpdateSource  (0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,E,Jp,Jm);
            
            
        }
        
        void Update(){
            Update(SUNcGaugeLinks::U,SUNcElectricFields::E,SUNcChargeDensity::Jp,SUNcChargeDensity::Jm);
        }
	   
    }
    
}


#endif
