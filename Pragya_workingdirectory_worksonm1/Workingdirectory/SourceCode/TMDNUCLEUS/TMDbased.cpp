#ifndef __TMD_BASED_DISTRIBUTIONS__
#define __TMD_BASED_DISTRIBUTIONS__

namespace TmdBasedCollision{
    
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
	
    void SetChargeDensities(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zStep,ChargeDensity *rho,ChargeDensity *CFs,char D,Nucleus *N){
		
	   INT z=zStep;
				
	   for(INT y=yLow;y<=yHigh;y++){
		  for(INT x=xLow;x<=xHigh;x++){
			 for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				 
				 double xInFm=(x-Lattice::N[0]/2)*Parameters::aXInFm;
				 double yInFm=(y-Lattice::N[1]/2)*Parameters::aYInFm;
				 double zInFm;
				 
				 if(D=='L'){
//					zInFm=(zStep-(Lattice::N[2]/4-1))*Parameters::aZInFm;
					zInFm=(zStep-785)*Parameters::aZInFm;
				 }

				 if(D=='R'){

					zInFm=(zStep-1261)*Parameters::aZInFm;
				 }
				 
				 // rho_lat= \sqrt(azInfm*T(z)/az)*Lattice::aCube //
				 rho->Get(x,y,z,a)[0]=Lattice::aCube*sqrt(Parameters::aZInFm*N->GetProfile(xInFm,yInFm,zInFm)/Lattice::a[2])*CFs->Get(x,y,z,a)[0];

				 FourierVariables->SetX(x,y,a,rho->GetValue(x,y,z,a));
			 }
		  }
	   }
		
	  //CUTTING OFF THE SMALL MODES BY TRANSFORMIMG \rho INTO \tilde{\rho}=\rho*k^2/k^2+m^2
      // PERFORM FFT X->P //
	  FourierVariables->ExecuteXtoP();
	  
	  DOUBLE pSqr; COMPLEX RhoTilde; DOUBLE NormalizationFactor=(1.0)/(Lattice::N[0]*Lattice::N[1]);
	  
	  for(INT pYIndex=0;pYIndex<=Lattice::N[1]-1;pYIndex++){
		 for(INT pXIndex=0;pXIndex<=Lattice::N[0]-1;pXIndex++){
			
			// DETERMINE MOMENTUM //
		    pSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pXIndex/DOUBLE(Lattice::N[0])))/SQR(Lattice::a[0])+(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pYIndex/DOUBLE(Lattice::N[1])))/SQR(Lattice::a[1]);

			for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
			   
			   // DETERMINE NEW RHO'S VALUE //
			   RhoTilde=NormalizationFactor*FourierVariables->GetP(pXIndex,pYIndex,a)*pSqr/(pSqr+SQR(Parameters::MEff));
			   
			   // SET VALUE //
			   FourierVariables->SetP(pXIndex,pYIndex,a,RhoTilde);
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
	 
		  
	  // SETTING RHO FOR SOLVING LAPLACE EQUATION  //
	  for(INT y=yLow;y<=yHigh;y++){
		 for(INT x=xLow;x<=xHigh;x++){
			for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
			   FourierVariables->SetX(x,y,a,rho->GetValue(x,y,z,a));
			}
		 }
	  }
   }

	   
    void SetChargeDensities(INT zStep,ChargeDensity *rho,ChargeDensity *CFs,char D,Nucleus *N){
	   
	   SetChargeDensities(0,Lattice::N[0]-1,0,Lattice::N[1]-1,zStep,rho,CFs,D,N);
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
				
				// COMPUTE WILSON LINE UPDATE //
				
				//////////////////////////////////////////////////////////////////////////
				// V = Product_{zLow,zHigh}exp(igA^-lat)  where A^-lattice =ga[2]A^-cont //
				//////////////////////////////////////////////////////////////////////////
				
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    ColorVectorBuffer[a]=real(FourierVariables->GetX(x,y,a));
				}
				 				
				if(D=='R'){
					SUNcAlgebra::Operations::MatrixIExp(+DOUBLE(1.0)/M_SQRT2,ColorVectorBuffer,WilsonLineUpdate);
				}
				
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
    
    
    void UpdateWilsonLines(INT zStep,ChargeDensity *rho,ChargeDensity *CFs,char D,WilsonLines *V,Nucleus *N){
	   
	   // SET COLOR CHARGE DENSITIES //
	   SetChargeDensities(zStep,rho,CFs,D,N);
	   
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
    
    
    void SetWilsonLines(ChargeDensity *rho,ChargeDensity *CFs,WilsonLines *V,char D){
	   
	   // COMMADNLINE OUTPUT //
	   std::cerr << "#SETTING WILSON LINES FOR "<<D<<" SITTING NUCLEUS" << std::endl;
	   
	   // SET TARGET WILSON LINE //
	   Wil=new SU_Nc_FUNDAMENTAL_FORMAT[SUNcGroup::MatrixSize*Lattice::N[0]*Lattice::N[1]];
	   
	   // SET WILSONS LINES TO UNIT MATRICES //
	   SetIdentity();
	   
	   // INITIALIZE FOURIER SPACE VARIABLES //
	   FourierVariables=new FFT2D(Lattice::N[0],Lattice::N[1],SUNcAlgebra::NumberOfGenerators);
		
	   Nucleus *N=new Nucleus(Parameters::AtomicNumber);
	
//	   N->OutputProfile(StringManipulation::StringCast(IO::OutputDirectory,"NucleonProfile",D,".txt").c_str(),Lattice::N[0],Lattice::N[1],Lattice::N[2],D);
//
//	   N->OutputNucleonPositions(StringManipulation::StringCast(IO::OutputDirectory,"NucleonPosition",D,".txt").c_str());
		
	   //RIGHT MEANS RIGHT POSITIONED
	   if(D=='R'){
		  
		  for(INT zStep=0;zStep<=Lattice::N[2]-1;zStep++){
			 
			 UpdateWilsonLines(zStep,rho,CFs,D,V,N);
		  }
		  
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;
		  
	   }
	   
	   //L MEANS LEFT POSITIONED
	   if(D=='L'){
		  
		  for(INT zStep=Lattice::N[2]-1;zStep>=0;zStep--){
			 
			 UpdateWilsonLines(zStep,rho,CFs,D,V,N);
		  }
		  // CLEAN UP FOURIER SPACE VARIABLES //
		  delete FourierVariables;
		  
	   }
		
		rho->SynchronizeGhostCells();
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

