#ifndef __GEN_COLLIDING_OBJECTS__
#define __GEN_COLLIDING_OBJECTS__

namespace CollidingObjects{
    
    SU_Nc_FUNDAMENTAL_FORMAT *Wil;
    
    FFT2D *FourierVariables;

    void SetIdentity(INT xLow,INT xHigh,INT yLow,INT yHigh){
	   
	   for(INT y=yLow;y<=yHigh;y++){
		  for(INT x=xLow;x<=xHigh;x++){
			 
			 COPY_SUNcMatrix(&Wil[SUNcMatrixIndex2D(x,y)],SUNcGroup::UnitMatrix);
		  }
	   }
    }
    
    void SetIdentity(){
	   
	   SetIdentity(0,Lattice::N[0]-1,0,Lattice::N[1]-1);
    }
    
	 /////////////////////////////////////////////////////////////////////////////////////////////////
	//   rho(z,x_perp)=g*a[2]*a[0]*a[1]*Qs0A*T(z)*ChargeFluctuation(x_perp)                 //
	//   T(z)==1/(sqrt(2*M_PI*Sigma^2) Exp(-1/2 *((z-z0)^2/Sigma^2)     //
    ////////////////////////////////////////////////////////////////////////////////////////////////

    DOUBLE LongitudinalProfile(INT z,INT Center,DOUBLE Width){
	   return (std::exp(-0.5*SQR((z-Center)/DOUBLE(Width))))/(sqrt(2*M_PI*Width*Width));
    }
	   
    void SetChargeDensities(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zStep,INT Center,DOUBLE* ChargeFluctuations,ChargeDensity *rhoOld,ChargeDensity *rho,char D){
	   
	   INT z=zStep; DOUBLE Value; DOUBLE Normalization=sqrt(Lattice::a[0]*Lattice::a[1]);
	   		
	   for(INT y=yLow;y<=yHigh;y++){
		  for(INT x=xLow;x<=xHigh;x++){
			 
			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				
				if(D=='L'){
			    	Value=Parameters::LeftNucleus*Normalization*LongitudinalProfile(zStep,Center,Parameters::SigmaLOverAzL)*(ChargeFluctuations[SUNcAlgebraIndex2D(x,y,a)]);
                }
                if(D=='R'){
				    Value=Parameters::RightNucleus*Normalization*LongitudinalProfile(zStep,Center,Parameters::SigmaLOverAzR)*(ChargeFluctuations[SUNcAlgebraIndex2D(x,y,a)]);
                }
				 
				rhoOld->Get(x,y,z,a)[0]=Value;
				FourierVariables->SetX(x,y,a,rhoOld->GetValue(x,y,z,a));
			 }
		  }
	   }
	   
	   // SMOOTHEN CHARGE DENSITY //
	   
	   // PERFORM FFT X->P //
	   FourierVariables->ExecuteXtoP();
	   
	   // SMOOTHEN IN FOURIER SPACE //
	   DOUBLE pSqr; COMPLEX SmoothValue; DOUBLE NormalizationFactor=(1.0)/(Lattice::N[0]*Lattice::N[1]);
	   
	   for(INT pYIndex=0;pYIndex<=Lattice::N[1]-1;pYIndex++){
		  for(INT pXIndex=0;pXIndex<=Lattice::N[0]-1;pXIndex++){
			 
			 // DETERMINE MOMENTUM //
			pSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pXIndex/DOUBLE(Lattice::N[0])))/SQR(Lattice::a[0])+(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pYIndex/DOUBLE(Lattice::N[1])))/SQR(Lattice::a[1]);

			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				
				// DETERMINE SMOOTH VALUE //
				SmoothValue=NormalizationFactor*FourierVariables->GetP(pXIndex,pYIndex,a)*(std::exp(-pSqr/(2*SQR(Parameters::Lambda))))*pSqr/(pSqr+SQR(Parameters::MEff));
				
				// SET VALUE //
				FourierVariables->SetP(pXIndex,pYIndex,a,SmoothValue);
			 }
		  }
	   }
	   
	   // PERFORM FFT P->X //
	   FourierVariables->ExecutePtoX();
	   
		  for(INT y=yLow;y<=yHigh;y++){
			 for(INT x=xLow;x<=xHigh;x++){
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

				    rho->Get(x,y,z,a)[0]=real(FourierVariables->GetX(x,y,a));
				}
			 }
		  }
	  
	   
	   // SETTING SMOOTHENED RHO FOR SOLVING LAPLACE EQUATION  //
	   for(INT y=yLow;y<=yHigh;y++){
		  for(INT x=xLow;x<=xHigh;x++){
			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				FourierVariables->SetX(x,y,a,rho->GetValue(x,y,z,a));
			 }
		  }
	   }	   
    }
    
    
    void SetChargeDensities(INT zStep,INT Center,DOUBLE* ChargeFluctuations,ChargeDensity *rhoOld,ChargeDensity *rho,char D){
	   
	   SetChargeDensities(0,Lattice::N[0]-1,0,Lattice::N[1]-1,zStep,Center,ChargeFluctuations,rhoOld,rho,D);
    }
    
    void SetChargeDensityByEvolvedRho(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zStep,ChargeDensity *J){
	   
	   INT z=zStep;
	   
	   for(INT y=yLow;y<=yHigh;y++){
		  for(INT x=xLow;x<=xHigh;x++){
			 
			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				
				FourierVariables->SetX(x,y,a,J->GetValue(x,y,zStep,a));
			 }
		  }
	   }
    }
    
    
    void SetChargeDensityByEvolvedRho(INT zStep,ChargeDensity *J){
	   
	   SetChargeDensityByEvolvedRho(0,Lattice::N[0]-1,0,Lattice::N[1]-1,zStep,J);
    }
    
    void SolveLaplaceEquation(INT pXLow,INT pXHigh,INT pYLow,INT pYHigh){
	   
#pragma omp parallel
	   {
		  DOUBLE pSqr; COMPLEX  LaplaceValue; DOUBLE NormalizationFactor=(1.0)/(Lattice::N[0]*Lattice::N[1]);
		  
#pragma omp for collapse(2)
		  for(INT pYIndex=pYLow;pYIndex<=pYHigh;pYIndex++){
			 for(INT pXIndex=pXLow;pXIndex<=pXHigh;pXIndex++){
				
				// DETERMINE MOMENTUM //
				pSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pXIndex/DOUBLE(Lattice::N[0])))/SQR(Lattice::a[0])+(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pYIndex/DOUBLE(Lattice::N[1])))/SQR(Lattice::a[1]);
				
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    
				    // DETERMINE LAPLACE //
				    LaplaceValue=NormalizationFactor*FourierVariables->GetP(pXIndex,pYIndex,a)/pSqr; 
				    
				    // EXCLUDE ZERO MODE //
				    if(pXIndex==0 && pYIndex==0){
					   FourierVariables->SetP(pXIndex,pYIndex,a,COMPLEX(0.0,0.0));
				    }
				    
				    // SET LAPLACE VALUE //
				    else{
					   FourierVariables->SetP(pXIndex,pYIndex,a,LaplaceValue);
				    }
				}
			 }
		  }
	   } // END PARALLEL
    }
    
    
    void SolveLaplaceEquation(){
	   
	   // PERFORM FFT X->P //
	   FourierVariables->ExecuteXtoP();
	   
	   // SOLVE LAPLACE EQUATION IN MOMENTUM SPACE //
	   SolveLaplaceEquation(0,Lattice::N[0]-1,0,Lattice::N[1]-1);
	   
	   // PERFORM FFT P->X //
	   FourierVariables->ExecutePtoX();
	   
    }
    
    
    void SetNewWilsonLines(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zStep,char D,WilsonLines *V){
	   
	   INT z=zStep;
	   
#pragma omp parallel
	   {
		  DOUBLE ColorVectorBuffer[SUNcAlgebra::NumberOfGenerators];
		  
		  SU_Nc_FUNDAMENTAL_FORMAT WilsonLineUpdate[SUNcGroup::MatrixSize];
		  SU_Nc_FUNDAMENTAL_FORMAT VOld[SUNcGroup::MatrixSize];
		  
#pragma omp for
		  for(INT y=yLow;y<=yHigh;y++){
			 for(INT x=xLow;x<=xHigh;x++){
				
				// COMPUTE WILSON LINE UPDATE //WilsonLines
				
				//////////////////////////////////////////////////////////////////////////
				// V = Product_{zLow,zHigh}exp(igA^-lat)  where A^-lattice =a[2]A^-cont //
				//////////////////////////////////////////////////////////////////////////

				
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    ColorVectorBuffer[a]=real(FourierVariables->GetX(x,y,a));
				}
				
				if(D=='R'){
				    SUNcAlgebra::Operations::MatrixIExp(+DOUBLE(1.0)/M_SQRT2,ColorVectorBuffer,WilsonLineUpdate);
				}
				//M_SQRT2
				if(D=='L'){
				    SUNcAlgebra::Operations::MatrixIExp(-DOUBLE(1.0)/M_SQRT2,ColorVectorBuffer,WilsonLineUpdate);
				}
				
				// COMPUTE NEW WILSON LINE //
				
				COPY_SUNcMatrix(VOld,&Wil[SUNcMatrixIndex2D(x,y)]);
				
				SUNcGroup::Operations::UU(WilsonLineUpdate,VOld,&Wil[SUNcMatrixIndex2D(x,y)]);
				
				COPY_SUNcMatrix(V->Get(x,y,zStep),&Wil[SUNcMatrixIndex2D(x,y)]);
				
			 }
		  }
		  
	   } // END PARALLEL
	   
    }
    
    
    void SetNewWilsonLines(INT zStep,char D,WilsonLines *V){
	   
	   SetNewWilsonLines(0,Lattice::N[0]-1,0,Lattice::N[1]-1,zStep,D,V);
    }
    
    
    void UpdateWilsonLines(INT zStep,INT Center, DOUBLE* ChargeFluctuations,ChargeDensity *rhoOld,ChargeDensity *rho,char D,WilsonLines *V){
	   
	   // SET COLOR CHARGE DENSITIES //
	   SetChargeDensities(zStep, Center, ChargeFluctuations,rhoOld,rho,D);
	   
	   // SOLVE LAPLACE EQUATION //
	   SolveLaplaceEquation();
	   
	   // UPDATE WILSON LINES //
	   SetNewWilsonLines(zStep,D,V);
	   
    }
    
    void UpdateWilsonLinesAtFirstTimeSlice(INT zStep,ChargeDensity *J,WilsonLines *V, char D){
	   
	   // SET COLOR CHARGE DENSITIES //
	   SetChargeDensityByEvolvedRho(zStep,J);
	   
	   // SOLVE LAPLACE EQUATION //
	   SolveLaplaceEquation();
	   
	   // UPDATE WILSON LINES //
	   SetNewWilsonLines(zStep,D,V);
	   
    }
    
	void SetWilsonLines(ChargeDensity *rhoOld,ChargeDensity *rho,WilsonLines *V,char D,INT Center){    
	   
	   // COMMADNLINE OUTPUT //
	   std::cerr << "#SETTING WILSON LINES FOR "<<D<<" SITTING NUCLEUS" << std::endl;
	   
	   // SET TARGET WILSON LINE //
	   Wil=new SU_Nc_FUNDAMENTAL_FORMAT[SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]];
	   
	   // SET WILSONS LINES TO UNIT MATRICES //
	   SetIdentity();
	   
	   // INITIALIZE FOURIER SPACE VARIABLES //
	   FourierVariables=new FFT2D(Lattice::N[0],Lattice::N[1],SUNcAlgebra::NumberOfGenerators);

	   DOUBLE *ChargeFluctuations;
	   ChargeFluctuations= new DOUBLE[SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]];
	   
	   for(INT y=0;y<=Lattice::N[1]-1;y++){
		  for(INT x=0;x<=Lattice::N[0]-1;x++){
			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

				ChargeFluctuations[SUNcAlgebraIndex2D(x,y,a)]=Parameters::Qs0A*RandomNumberGenerator::Gauss();
			 }
		  }
	   }
	   
	   //RIGHT MEANS RIGHT POSITIONED
	   if(D=='R'){
		  
		  for(INT zStep=0;zStep<=Lattice::N[2]-1;zStep++){
			 
			 UpdateWilsonLines(zStep,Center,ChargeFluctuations,rhoOld,rho,D,V);
		  }
		  
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;

		  delete [] ChargeFluctuations;
	   }
	   
	   //L MEANS LEFT POSITIONED
	   if(D=='L'){
		  
		  for(INT zStep=Lattice::N[2]-1;zStep>=0;zStep--){
			 
			 UpdateWilsonLines(zStep,Center,ChargeFluctuations,rhoOld,rho,D,V);
		  }
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;

		  delete [] ChargeFluctuations;
	   }
	   
    }



	void SetWilsonLines(ChargeDensity *rhoOld,ChargeDensity *rho,WilsonLines *V,char D,INT Center, DOUBLE QsA){
	   
	   // COMMADNLINE OUTPUT //
	   std::cerr << "#SETTING WILSON LINES FOR "<<D<<" SITTING NUCLEUS" << std::endl;
	   
	   // SET TARGET WILSON LINE //
	   Wil=new SU_Nc_FUNDAMENTAL_FORMAT[SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]];
	   
	   // SET WILSONS LINES TO UNIT MATRICES //
	   SetIdentity();
	   
	   // INITIALIZE FOURIER SPACE VARIABLES //
	   FourierVariables=new FFT2D(Lattice::N[0],Lattice::N[1],SUNcAlgebra::NumberOfGenerators);

	   DOUBLE *ChargeFluctuations;
	   ChargeFluctuations= new DOUBLE[SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]];
	   
	   for(INT y=0;y<=Lattice::N[1]-1;y++){
		  for(INT x=0;x<=Lattice::N[0]-1;x++){
			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){

				ChargeFluctuations[SUNcAlgebraIndex2D(x,y,a)] = QsA * RandomNumberGenerator::Gauss();
			 }
		  }
	   }
	   
	   //RIGHT MEANS RIGHT POSITIONED
	   if(D=='R'){
		  
		  for(INT zStep=0;zStep<=Lattice::N[2]-1;zStep++){
			 
			 UpdateWilsonLines(zStep,Center,ChargeFluctuations,rhoOld,rho,D,V);
		  }
		  
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;

		  delete [] ChargeFluctuations;
	   }
	   
	   //L MEANS LEFT POSITIONED
	   if(D=='L'){
		  
		  for(INT zStep=Lattice::N[2]-1;zStep>=0;zStep--){
			 
			UpdateWilsonLines(zStep,Center,ChargeFluctuations,rhoOld,rho,D,V);
		  }
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;

		  delete [] ChargeFluctuations;
	   }
	   
    }




    
    void SetWilsonLinesForFirstStep(WilsonLines *V,ChargeDensity *J,char D){
	   
	   // COMMADNLINE OUTPUT //
	   std::cerr << "#SETTING WILSON LINES FOR "<< D<<" SITING NUCLEUS AGAIN FOR ELECTRIC FIELDS "<<std::endl;
	   
	   // SET TARGET WILSON LINE //
	   Wil=new SU_Nc_FUNDAMENTAL_FORMAT[SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]];
	   
	   // SET WILSONS LINES TO UNIT MATRICES //
	   SetIdentity();
	   
	   // INITIALIZE FOURIER SPACE VARIABLES //
	   FourierVariables=new FFT2D(Lattice::N[0],Lattice::N[1],SUNcAlgebra::NumberOfGenerators);
	   
	   //R MEANS RIGHT POSITIONED
	   if(D=='R'){
		  
		  for(INT zStep=0;zStep<=Lattice::N[2]-1;zStep++){
			 UpdateWilsonLinesAtFirstTimeSlice(zStep,J,V,D);
		  }
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;
	   }
	   
	   //L MEANS LEFT POSITIONED
	   if(D=='L'){
		  
		  for(INT zStep=Lattice::N[2]-1;zStep>=0;zStep--){
			 
			 UpdateWilsonLinesAtFirstTimeSlice(zStep,J,V,D);
		  }
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;
	   }
    }
    
}

#endif

