namespace GaussLawRestoration{
    
    //GAUSS VIOLATION INDEX -- MACRO TO RETURN THE POSITION OF A CERTAIN COMPONENT OF COLORED FIELD
#define GaussViolationIndex(x,y,z,a) ((a)+SUNcAlgebra::NumberOfGenerators*(E->Index3D((x),(y),(z))))
    
    // GAUSS VIOLATION ARRAY //
    DOUBLE *GaussViolation;
    
    //GLOBAL MAXIMUM AND AVERAGE VIOLATION OF GAUSS LAW
    DOUBLE GlobalMaxViolation;
    DOUBLE GlobalAvgViolation;
    
    
    void UpdateDeviation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp, ChargeDensity *Jm){
	   
	   
	   //MAXIMUM VIOLATION
	   DOUBLE MaxViolation=0.0;
	   
	   //SQR SUM OF ALL VIOLATIONS
	   DOUBLE SqrSumViolation=0.0;
	   
	   
#pragma omp parallel
	   {
		  
		  //ALLOCATE BUFFERS
		  SET_GAUSS_LAW_BUFFERS();
		  
		  //LOCAL VIOLATION
		  DOUBLE LocalViolation[SUNcAlgebra::NumberOfGenerators];
		  
		  // UPDATE ELECTRIC FIELDS //
#pragma omp for reduction( + : SqrSumViolation)  reduction ( max : MaxViolation)
		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
				    
				    // GET LOCAL VIOLATION //
				    COMPUTE_GAUSS_VIOLATION(LocalViolation,x,y,z);
				    
				    // COMPUTEVALL COMPONENTS //
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   
					   GaussViolation[GaussViolationIndex(x,y,z,a)]=LocalViolation[a];
					   
					   //UPDATE MAXIMUM
					   MaxViolation=std::max(MaxViolation,std::abs(GaussViolation[GaussViolationIndex(x,y,z,a)]));
					   
					   //UPDATE SQR SUM
					   SqrSumViolation+=SQR(GaussViolation[GaussViolationIndex(x,y,z,a)]);
				    }
				    
				}
			 }
		  }
		  
	   } // END PARALLEL
	   
	   //GET GLOBAL RESULTS
	   GlobalMaxViolation=MaxViolation;
	   GlobalAvgViolation=sqrt(SqrSumViolation/(SUNcAlgebra::NumberOfGenerators*U->Volume));
    
    }

    void UpdateDeviation(GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp,ChargeDensity *Jm){
	   
	   std::cerr << "# NOT SURE HERE" << std::endl;
	   UpdateDeviation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E,Jp,Jm);
	   
    }
    
    
    // UPDATE ELECTRIC FIELDS //
    void UpdateElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,DOUBLE Gamma){
	   
#pragma omp parallel
	   {
		  
		  //ALLOCATE BUFFERS
		  SU_Nc_ALGEBRA_FORMAT GaussXUp[SUNcAlgebra::NumberOfGenerators];
		  SU_Nc_ALGEBRA_FORMAT GaussYUp[SUNcAlgebra::NumberOfGenerators];
		  SU_Nc_ALGEBRA_FORMAT GaussZUp[SUNcAlgebra::NumberOfGenerators];
		  
		  // UPDATE ELECTRIC FIELDS //
#pragma omp for
		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
				    
				    // GET TRANSPORTED UPPER VALUES //
				    SUNcAlgebra::Operations::AdjointMultiplication(U->Get(x,y,z,0),&GaussViolation[GaussViolationIndex(x+1,y,z,0)],GaussXUp);
				    SUNcAlgebra::Operations::AdjointMultiplication(U->Get(x,y,z,1),&GaussViolation[GaussViolationIndex(x,y+1,z,0)],GaussYUp);
				    SUNcAlgebra::Operations::AdjointMultiplication(U->Get(x,y,z,2),&GaussViolation[GaussViolationIndex(x,y,z+1,0)],GaussZUp);
				    
				    // UPDATE ELECTRIC FIELDS //
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   
					   E->Get(x,y,z,0,a)[0]-=Gamma*(GaussViolation[GaussViolationIndex(x,y,z,a)]-GaussXUp[a]);
					   E->Get(x,y,z,1,a)[0]-=Gamma*(GaussViolation[GaussViolationIndex(x,y,z,a)]-GaussYUp[a]);
					   E->Get(x,y,z,2,a)[0]-=Gamma*(GaussViolation[GaussViolationIndex(x,y,z,a)]-GaussZUp[a]);
					   
				    }
				    
				    
				}
				
			 }
		  }
		  
		  
	   } // END PARALLEL
	   
	   // SYNCHRONIZE BOUNDARIES //
    }
    
    void UpdateElectricFields(GaugeLinks *U,ElectricFields *E,DOUBLE Gamma){
	   
	   UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E,Gamma);

	   
    }
    
    
    // COMPLETE UPDATE STEP //
    void UpdateStep(GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp, ChargeDensity *Jm,DOUBLE Gamma){
	   
	   UpdateDeviation(U,E,Jp,Jm);
	   UpdateElectricFields(U,E,Gamma);
	   
    }
    
    
    void Restore(GaugeLinks *U,ElectricFields *E,ChargeDensity *Jp, ChargeDensity *Jm){
	   
	   // ALLOCATE MEMORY //
	   GaussViolation=new DOUBLE[U->Volume*SUNcAlgebra::NumberOfGenerators];
	   
	   // STEP COUNTER //
	   INT nSteps=0;
	   
	   //RESET BREAK CONDITION //
	   GlobalMaxViolation=1;
	   
	   std::cerr << "#GAUSS LAW FIXING -- STEPS MAX-DEV AVG-DEV" << std::endl;
	   
	   // CHECK INITIAL DEVIATION //
	   UpdateDeviation(U,E,Jp,Jm);
	   std::cerr << nSteps << " " << GlobalMaxViolation << " " << GlobalAvgViolation << std::endl;
	   
	   // ITERATIVE UPDATE //
	   while(nSteps<30000 && GlobalMaxViolation>std::pow(10.0,-MAX_DIGITS_PRECISION+1)){
		  
		  UpdateStep(U,E,Jp,Jm,7.0/24.0);
		  UpdateStep(U,E,Jp,Jm,7.0/48.0);
		  UpdateStep(U,E,Jp,Jm,7.0/72.0);

		  nSteps++;
		  
		  if(nSteps%100==0){
			 
			 std::cerr << nSteps << " " << GlobalMaxViolation << " " << GlobalAvgViolation << std::endl;
			 
		  }
		  
		  
	   }
	   
	   std::cerr << "#GAUSS LAW COMPLETED  -- STEPS MAX-DEV AVG-DEV" << std::endl;
	   
	   // CHECK INITIAL DEVIATION //
	   UpdateDeviation(U,E,Jp,Jm);
	   
	   std::cerr << nSteps << " " << GlobalMaxViolation << " " << GlobalAvgViolation << std::endl;
	   
	   // FREE MEMORY //
	   delete[] GaussViolation;
    }
    
    void Restore(){
	   Restore(SUNcGaugeLinks::U,SUNcElectricFields::E,SUNcChargeDensity::JpStat,SUNcChargeDensity::JmStat);
    }
    
    
}
