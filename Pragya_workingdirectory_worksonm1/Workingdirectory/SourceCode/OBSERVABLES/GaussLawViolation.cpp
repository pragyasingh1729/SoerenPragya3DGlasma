#ifndef __GAUSSLAW__CPP__
#define __GAUSSLAW__CPP__

namespace Observables {
  
    namespace GaussLaw{

	   //GLOBAL MAXIMUM AND AVERAGE VIOLATION OF GAUSS LAW
	   
	   DOUBLE GlobalMaxViolation;
	   
	   DOUBLE GlobalAvgViolation;
	   
	   //COMPUTES THE VIOLATION OF THE GAUSS LAW CONSTRAIN
	   
	   void UpdateViolation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,ChargeDensity *JpStat, ChargeDensity *JmStat,std::string Gfname){
		  
		  // OUTPUT FILES //
		  std::string GOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Gfname,"ID",RandomNumberGenerator::MySEED,".txt");
		  
		  // OPEN FILES //
		  std::ofstream GOutStream;
		  GOutStream.open(GOutputFile.c_str());
		  
		  // SET PRECISION //
		  GOutStream.precision(OUTPUT_PRECISION);
		  
		  //ALLOCATE BUFFERS
		  SET_GAUSS_LAW_BUFFERS();
		  
		  //LOCAL VIOLATION
		  DOUBLE LocalViolation[SUNcAlgebra::NumberOfGenerators];
		   
		   SUNcChargeDensity::JpStat->SynchronizeGhostCells();
		   SUNcChargeDensity::JmStat->SynchronizeGhostCells();
		   
		   SUNcGaugeLinks::U->SynchronizeGhostCells();
		   SUNcElectricFields::E->SynchronizeGhostCells();
	  
		  //COMPUTE LOCAL VIOLATION AT EACH POINT
		  for(INT z=zLow;z<=zHigh;z++){
			 
			 DOUBLE AvgGaussOverTransProfile=0.0;
			 
			 for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
				    
				    // GET LOCAL VIOLATION //
				    COMPUTE_GAUSS_VIOLATION(LocalViolation,x,y,z);
				    
				    // ERROR OUTPUT //
				    double Violation=0.0;
				    
				    for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   
					   Violation+=SQR(LocalViolation[a]);
				    }
				    
				    AvgGaussOverTransProfile+=Violation;
				}
			 }
			 
			 GOutStream << z << "  " << AvgGaussOverTransProfile/(Lattice::N[0]*Lattice::N[1]) << std::endl;
			 
			 //GlobalMaxViolation=MaxViolation;
			 //GlobalAvgViolation=sqrt(SqrSumViolation/(SUNcAlgebra::NumberOfGenerators*U->Volume));
		  }
		  
		  
		  // CLOSE OUTPUT STREAM //
		  GOutStream.close();
	   }
	   	   
	   
	   void UpdateViolation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,ChargeDensity *JpStat, ChargeDensity *JmStat){
		  
	       //MAXIMUM VIOLATION
		  DOUBLE MaxViolation=0.0;
		  
		  //SQR SUM OF ALL VIOLATIONS//
		  DOUBLE SqrSumViolation=0.0;
		  
#pragma omp parallel
		  
		  {
			 //ALLOCATE BUFFERS
			 
			 SET_GAUSS_LAW_BUFFERS();
			 //LOCAL VIOLATION
			 
			 DOUBLE LocalViolation[SUNcAlgebra::NumberOfGenerators];
			 
			 //COMPUTE LOCAL VIOLATION AT EACH POINT
			 
#pragma omp for reduction( + : SqrSumViolation) reduction( max : MaxViolation)
			 
			 for(INT z=zLow;z<=zHigh;z++){
				for(INT y=yLow;y<=yHigh;y++){
				    for(INT x=xLow;x<=xHigh;x++){
					  
					   // GET LOCAL VIOLATION //
					   COMPUTE_GAUSS_VIOLATION(LocalViolation,x,y,z);
					   
					   //COMPUTE VIOLATION FOR EACH COLOR COMPONENT
					   for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						  
					     //UPDATE MAXIMUM
						  MaxViolation=std::max(MaxViolation,std::abs(LocalViolation[a]));
						  //UPDATE SQR SUM
						  
						  SqrSumViolation+=SQR(LocalViolation[a]);
					   }
				    }
				}
			 }
		  } // END PARALLEL
		  
		  //GET GLOBAL RESULTS
		  
		  GlobalMaxViolation=MaxViolation;
		  
		  GlobalAvgViolation=sqrt(SqrSumViolation/(SUNcAlgebra::NumberOfGenerators*U->Volume));
		  
	   }
	   
	   void CheckViolation(GaugeLinks *U,ElectricFields *E,ChargeDensity *JpStat, ChargeDensity *JmStat,std::string Gfname){
		  
		  // CHECK VIOLATION OF GAUSS LAW //
		  UpdateViolation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E,JpStat,JmStat,Gfname);
		  
		  // COMMANDLINE OUTPUT //
		  std::cerr << "#GAUSS LAW VIOLATION AT T=" << Dynamics::Time() << " MAX="  << GlobalMaxViolation << " AVG=" << GlobalAvgViolation << std::endl;

	   }
	   
		void CheckViolation(std::string Gfname){
		   
		   CheckViolation(SUNcGaugeLinks::U,SUNcElectricFields::E,SUNcChargeDensity::JpStat,SUNcChargeDensity::JmStat,Gfname);
		}
	
	   void SaveViolation(std::string GVfname,GaugeLinks *U,ElectricFields *E,ChargeDensity *JpStat, ChargeDensity *JmStat){
		  	
		  std::ofstream GVOutStream;
		  // OUTPUT FILES //
		  
		  std::string GVOutputFile=StringManipulation::StringCast(IO::OutputDirectory,GVfname,RandomNumberGenerator::MySEED,".txt");
		  
		  // OPEN FILES //
		  GVOutStream.open(GVOutputFile.c_str(),std::ios_base::app);
		  
		  // SET PRECISION //
		  
		  GVOutStream.precision(OUTPUT_PRECISION);
		  
		  // CHECK VIOLATION OF GAUSS LAW //
		  
		  UpdateViolation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E,JpStat,JmStat);
		  
		  // COMMANDLINE OUTPUT //
		  
		  GVOutStream << Dynamics::Time() << " " << GlobalMaxViolation << "  " << GlobalAvgViolation << std::endl;
		  
		  GVOutStream.close();
	   }
	   
	   void SaveViolation(std::string GVfname){
		  
		  SaveViolation(GVfname,SUNcGaugeLinks::U,SUNcElectricFields::E,SUNcChargeDensity::JpStat,SUNcChargeDensity::JmStat);
	   }
    }
}

#endif


