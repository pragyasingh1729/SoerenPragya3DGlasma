#ifndef __SUNC_DYNAMICS_CPP__
#define __SUNC_DYNAMICS_CPP__


namespace Dynamics{
    
    //////////////
    //   TIME   //
    //////////////

	char InputOne[256]="";

	char InputTwo[256]="";

	char InputThree[256]="";


	DOUBLE LoadTime=0.0;
	
    FFT1D *FourierVariables;

    DOUBLE t0=0.0;
    
    //DISCRETIZED TIME AND NUMBER OF TIME STEPS
    DOUBLE t;
    
    INT tSteps=0;
    
    //TIME INCREMENT
	DOUBLE dTStep=0.08; //0.00044;
    
    //GET EVOLUTION TIME
    DOUBLE Time(){
	   return t;
    }
    
    
    ////////////////
    //   INIT    //
    ////////////////
    
	void LoadedTime(double T){
		tSteps=t/dTStep;  t=T;
		
	}

    void Reset(){
	   
	   t=t0;   tSteps=0;
    }
    
    namespace SUNc{
	   
	   //////////////////////////////////////////////////////////////////////////////////////////////////
	   //                         COMPUTE UPDATE OF THE LATTICE GAUGE LINKS                            //
	   //////////////////////////////////////////////////////////////////////////////////////////////////
	   //                                                                                              //
	   //   U_{mu}(x+dTau)= Exp(-ig a_{mu}^2/a^3 -g_{\mu\nu} E^{\nu} dTau / sqrt(-g)) U_{mu}(x)          //
	   //                                                                                              //
	   //////////////////////////////////////////////////////////////////////////////////////////////////
	   
	   void UpdateGaugeLinks(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,DOUBLE dTStep){
		  
		  //SET CONSTANTS CONSTANTS//
		  DOUBLE cU[Lattice::Dimension];
		  
		  for(INT mu=0;mu<Lattice::Dimension;mu++){
			 cU[mu]=-dTStep*SQR(U->a[mu])/(U->aCube);
		  }
		  
#pragma omp parallel
		  
		  {
			 //MATRIX BUFFERS
			 SU_Nc_FUNDAMENTAL_FORMAT LinkUpdate[SUNcGroup::MatrixSize];
			 SU_Nc_FUNDAMENTAL_FORMAT OldLink[SUNcGroup::MatrixSize];
			 
			 //UPDATE ALL GAUGE LINKS
#pragma omp for
			 
			 for(INT z=zLow;z<=zHigh;z++){
				for(INT y=yLow;y<=yHigh;y++){
				    for(INT x=xLow;x<=xHigh;x++){
					   
					   //UPDATE ALL LORENTZ COMPONENTS
					   for(INT mu=0;mu<Lattice::Dimension;mu++){
						  
						  //COMPUTE MATRIX EXPONENTIAL
						  SUNcAlgebra::Operations::MatrixIExp(cU[mu],E->Get(x,y,z,mu,0),LinkUpdate);
						  
						  //COPY THE OLD LINK
						  COPY_SUNcMatrix(OldLink,U->Get(x,y,z,mu));
						  
						  //COMPUTE UPDATED LINK AT TIME t+dt
						  SUNcGroup::Operations::UU(LinkUpdate,OldLink,U->Get(x,y,z,mu));
						   
					   }
					}
				}
			 }
			 
		  } // END PARALLEL
		  
		  U->SynchronizeGhostCells();
		  
	   }
		
		void SaveOldLink(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
			
		    SUNcGaugeLinks::U->SynchronizeGhostCells();
		
		#pragma omp parallel for schedule(static) collapse(3)
			for(INT z=zLow;z<=zHigh;z++){
				for(INT y=yLow;y<=yHigh;y++){
					 for(INT x=xLow;x<=xHigh;x++){
						 
						//UPDATE ALL LORENTZ COMPONENTS
						for(INT mu=0;mu<Lattice::Dimension;mu++){
							
							//SAVE THE OLD LINK AT TIME t to GET MAGNETIC FIELD AT B(t)
							COPY_SUNcMatrix(SUNcGaugeLinks::UOld->Get(x,y,z,mu),SUNcGaugeLinks::U->Get(x,y,z,mu));
						}
					}
				}
			}

			SUNcGaugeLinks::UOld->SynchronizeGhostCells();
		}
	   
	   //////////////////////////////////////////////////////////
	   //      COMPUTE UPDATE OF THE LATTICE GAUGE LINKS       //
	   //////////////////////////////////////////////////////////
	   
	   void UpdateElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp, ChargeDensity *Jm,DOUBLE dTStep){
		  
		  ///////////////////////////////////
		  // CLASSICAL YANG-MILLS DYNAMICS //
		  ///////////////////////////////////
		  
		  // SET DYNAMICAL CONSTANTS -2d\tau \sqrt{-g} a^{3} g^{\mu\alpha} g^{\nu\alpha}/(a_{\mu}^2 a_{\nu}^2) FOR GAUGE FORCE //
		  
		  DOUBLE cE[Lattice::Dimension];
		  
		  DOUBLE gamma=-2.0*dTStep*U->aCube;
		  
		  for(INT mu=0;mu<Lattice::Dimension;mu++){
			 
			 cE[mu]=1.0/SQR(E->a[mu]);
		  }
		  
		  SUNcGaugeLinks::U->SynchronizeGhostCells();
#pragma omp parallel
		  {
			 
			 //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
			 SET_ELEMENTARY_PLAQUETTE_BUFFERS();
			 SET_NEIGHBORING_PLAQUETTE_BUFFERS();
			 
			 // BUFFERS TO COMPUTE TRACES OF PLAQUETTES
			 SET_ELEMENTARY_COLOR_TRACE_BUFFERS();
			 SET_NEIGHBORING_COLOR_TRACE_BUFFERS();
			 
			 //UPDATE ALL ELECTRIC FIELDS  LINKS
#pragma omp for
			 
			 for(INT z=zLow;z<=zHigh;z++){
				for(INT y=yLow;y<=yHigh;y++){
				    for(INT x=xLow;x<=xHigh;x++){
					   
					   //COMPUTE ELEMENTARY PLAQUETTES AND COLOR TRACES
					   COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z,U);
					   COMPUTE_ELEMENTARY_COLOR_TRACES();
					   
					   //COMPUTE NEIGHBORING PLAQUETTES AND COLOR TRACES
					   COMPUTE_NEIGHBORING_PLAQUETTES(x,y,z,U);
					   COMPUTE_NEIGHBORING_COLOR_TRACES();
					   
					   //UPDATE ELECTRIC FIELD VARIABLES
					   for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						  
						  E->Get(x,y,z,0,a)[0]+=gamma*cE[0]*(cE[1]*(ReTrITaUxy[a]-ReTrITaUxMy[a])-cE[2]*(ReTrITaUzx[a]-ReTrITaUMzx[a]));
						  E->Get(x,y,z,1,a)[0]+=gamma*cE[1]*(cE[2]*(ReTrITaUyz[a]-ReTrITaUyMz[a])-cE[0]*(ReTrITaUxy[a]-ReTrITaUMxy[a]));
						  E->Get(x,y,z,2,a)[0]+=gamma*cE[2]*(cE[0]*(ReTrITaUzx[a]-ReTrITaUzMx[a])-cE[1]*(ReTrITaUyz[a]-ReTrITaUMyz[a]));

					   }

					   // CALCULATE CHANGE DUE TO CURRENT
					   for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

						  DOUBLE Jz=(Jp->GetValue(x,y,z,a)-Jm->GetValue(x,y,z,a))/M_SQRT2;

						  E->Get(x,y,z,2,a)[0]+=-(dTStep/Lattice::a[2])*Jz;
					   }

				    }
				}
			 }
			 
		  } // END PARALLEL
		  
		  E->SynchronizeGhostCells();
		  
	   }
		
	    void UpdateJStat(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ChargeDensity *JpNew,ChargeDensity *JmNew,ChargeDensity *JpOld,ChargeDensity *JmOld,DOUBLE dTStep){

		   #pragma omp parallel
		   {

			  SU_Nc_ALGEBRA_FORMAT JpDown[SUNcAlgebra::NumberOfGenerators];
			  SU_Nc_ALGEBRA_FORMAT JmDown[SUNcAlgebra::NumberOfGenerators];


			  #pragma omp for
				 for(INT z=zLow;z<=zHigh;z++){
					for(INT y=yLow;y<=yHigh;y++){
					    for(INT x=xLow;x<=xHigh;x++){

						   SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x,y,z-1,2),JpOld->Get(x,y,z-1,0),JpDown);
						   SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x,y,z-1,2),JmOld->Get(x,y,z-1,0),JmDown);

						   for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

							  JpNew->Get(x,y,z,a)[0]+= -(dTStep/Lattice::a[2])*(JpOld->GetValue(x,y,z,a)-JpDown[a]);
							  JmNew->Get(x,y,z,a)[0]+=  (dTStep/Lattice::a[2])*(JmOld->GetValue(x,y,z,a)-JmDown[a]);

						 }
					 }
				  }
			   }

		   }
		   
		   JmNew->SynchronizeGhostCells();
		   JpNew->SynchronizeGhostCells();
	    }
	   
	    void UpdateJ(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ChargeDensity *JpNew,ChargeDensity *JmNew, ChargeDensity *JpOld,ChargeDensity *JmOld,DOUBLE dTStep){

		   #pragma omp parallel
		   {

			  SU_Nc_ALGEBRA_FORMAT JpUp[SUNcAlgebra::NumberOfGenerators];
			  SU_Nc_ALGEBRA_FORMAT JmUp[SUNcAlgebra::NumberOfGenerators];


			  #pragma omp for
			  for(INT z=zLow;z<=zHigh;z++){
				 for(INT y=yLow;y<=yHigh;y++){
					for(INT x=xLow;x<=xHigh;x++){

					    SUNcAlgebra::Operations::AdjointMultiplication(U->Get(x,y,z,2),JpOld->Get(x,y,z+1,0),JpUp);
					    SUNcAlgebra::Operations::AdjointMultiplication(U->Get(x,y,z,2),JmOld->Get(x,y,z+1,0),JmUp);

					    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

						   JpNew->Get(x,y,z,a)[0]+= -(dTStep/Lattice::a[2])*(JpUp[a]-JpOld->GetValue(x,y,z,a));
						   JmNew->Get(x,y,z,a)[0]+=  (dTStep/Lattice::a[2])*(JmUp[a]-JmOld->GetValue(x,y,z,a));
					    }
					}
				 }
			  }

		   }
		   
		   JmNew->SynchronizeGhostCells();
		   JpNew->SynchronizeGhostCells();

	    }
	   
			 ////////////////////////////////////////////////////////////
			 // EVOLVING CURRENT VIA FOURIER TRANSFORMATION           //
			 ////////////////////////////////////////////////////////////

	   void EvolveRhoViaFourier(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ChargeDensity *JNew,ChargeDensity *JOld,DOUBLE dTStep,char D){

		  FourierVariables=new FFT1D(Lattice::N[2],SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]);
		  DOUBLE NormalizationFactor=(1.0)/(Lattice::N[2]);

		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

					   FourierVariables->SetX(z,SUNcAlgebraIndex2D(x,y,a),JOld->GetValue(x,y,z,a));
				    }
				}
			 }
		  }

		  FourierVariables->ExecuteXtoP();

		  
		  for(INT pZIndex=0;pZIndex<=Lattice::N[2]-1;pZIndex++){

			 DOUBLE OmegaSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2])))/SQR(Lattice::a[2]);

			 for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

					   COMPLEX JOldTilde=NormalizationFactor*FourierVariables->GetP(pZIndex,SUNcAlgebraIndex2D(x,y,a));

					   COMPLEX JNewTilde; 
					   
					   //L MEANS L POSITIONED
					   if(D=='L'){
						  
						  if(pZIndex<Lattice::N[2]/2){
							 JNewTilde = (std::exp(-ComplexI*sqrt(OmegaSqr)*dTStep)*JOldTilde);
						  }
						  else{
							 JNewTilde = (std::exp(+ComplexI*sqrt(OmegaSqr)*dTStep)*JOldTilde);
						  }
					   }
					   
					   //R MEANS R POSITIONED
					   if(D=='R'){
						  
						  if(pZIndex<Lattice::N[2]/2){
							 JNewTilde = (std::exp(+ComplexI*sqrt(OmegaSqr)*dTStep)*JOldTilde);
						  }
						  else{
							 JNewTilde = (std::exp(-ComplexI*sqrt(OmegaSqr)*dTStep)*JOldTilde);
						  }
					   }
					   
					   FourierVariables->SetP(pZIndex,SUNcAlgebraIndex2D(x,y,a),JNewTilde);
				    }
				}
			 }
		  }

		  FourierVariables->ExecutePtoX();

		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

				        JNew->Get(x,y,z,a)[0]=real(FourierVariables->GetX(z,SUNcAlgebraIndex2D(x,y,a)));
				    }
				}
			 }
		  }

		  delete FourierVariables;
		  JNew->SynchronizeGhostCells();

	   }
	   
	   // GLOBAL UPDATE INCLUDING SU(Nc) ONLY //
	   void Update(GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp, ChargeDensity *Jm, ChargeDensity *JpStat, ChargeDensity *JmStat){
		  
		  // EVOLVE ELECTRIC FIELD AT HALF-INTEGER TIME STEP -- E(t-dt/2)-> E(t+dt/2) //
		  Dynamics::SUNc::UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E,Jp,Jm,dTStep);

		  // EVOLVE CURRENTS J(t-dt/2) -> J(t+dt/2) //
		  Dynamics::SUNc::UpdateJStat(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,U,JpStat,JmStat,Jp,Jm,dTStep);
//		  Dynamics::SUNc::UpdateJStat(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,U,JmStat,Jm,dTStep);

		  // EVOLVE GAUGE LINKS FOR HALF STEP -- U(t)->U(t+dt/2) //
		  Dynamics::SUNc::UpdateGaugeLinks(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E,DOUBLE(0.5)*dTStep);
		   
		   //CACULATES MAGNETIC FIELDS AT t+dt/2 WHERE ELECTRIC FIELD EXISTS //
		  SaveOldLink(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1);
		  
		  // EVOLVE CURRENTS J(t) -> J(t+dt) //
		  Dynamics::SUNc::UpdateJ(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,U,Jp,Jm,JpStat,JmStat,dTStep);
//		  Dynamics::SUNc::UpdateJ(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,U,Jm,JmStat,dTStep);
		  
		  // EVOLVE GAUGE LINKS FOR HALF STEP -- U(t+dt/2)->U(t+dt) //
		  Dynamics::SUNc::UpdateGaugeLinks(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E,DOUBLE(0.5)*dTStep);

		  // INCREASE STEP COUNTER //
		  tSteps++;   t+=dTStep;
		  
	   }
	   
    }
    
}
#endif


