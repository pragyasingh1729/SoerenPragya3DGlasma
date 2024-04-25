#ifndef __SAVECONFIGURATION__CPP__
#define __SAVECONFIGURATION__CPP__

namespace IO{
    
    void SaveConfiguration(std::string Ufname,std::string Efname,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
	   // OUTPUT STREAMS //
	   std::ofstream UOutStream,EOutStream;
	   
	   // OUTPUT FILES //
	   std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Ufname,"ID",RandomNumberGenerator::MySEED,".txt");
	   std::string EOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Efname,"ID",RandomNumberGenerator::MySEED,".txt");
	   
	   // OPEN FILES //
	   UOutStream.open(UOutputFile.c_str());
	   EOutStream.open(EOutputFile.c_str());
	   
	   // SET PRECISION //
	   UOutStream.precision(OUTPUT_PRECISION);
	   EOutStream.precision(OUTPUT_PRECISION);
	   
	   // CREATE OUTPUT //
	   for(INT z=zLow;z<=zHigh;z++){
		  for(INT y=yLow;y<=yHigh;y++){
			 for(INT x=xLow;x<=xHigh;x++){
				for(INT mu=0;mu<Lattice::Dimension;mu++){
				    
				    // OUTPUT GAUGE LINKS //
				    UOutStream << x << " " << y << " " << z << " " << mu << " " << SUNcGroup::IO::MatrixToString(SUNcGaugeLinks::U->Get(x,y,z,mu))<<std::endl;
				    
				    // OUTPUT ELECTRIC FIELDS //
				    EOutStream << x << " " << y << " " << z<<" "<< mu;
				    
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   EOutStream << " " << SUNcElectricFields::E->Get(x,y,z,mu,a)[0];
				    }
				    
				    EOutStream << std::endl;
				    
				}
			 }
		  }
	   }
	   
	   // CLOSE OUTPUT STREAM //
	   UOutStream.close();
	   EOutStream.close();
	   
    }
    
    
    void SaveConfiguration(std::string Ufname,std::string Efname){
	   
	   #if(BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)
	   SaveConfiguration(Ufname,Efname,0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,-1,SUNcGaugeLinks::U->N[2]);
	   #endif
	   
	   #if(BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   SaveConfiguration(Ufname,Efname,0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1);
	   #endif
	   
    }
    
    
    void SaveConfigurationJ(std::string Jmfname){
	   
	   std::ofstream JmOutStream;
	   
	   // OUTPUT FILES //
	   std::string JmOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Jmfname,"ID",RandomNumberGenerator::MySEED,".txt");
	   
	   // OPEN FILES //
	   JmOutStream.open(JmOutputFile.c_str());
	   
	   // SET PRECISION //
	   JmOutStream.precision(OUTPUT_PRECISION);
	   
	   // CREATE OUTPUT //
	   for(INT z=-1;z<=SUNcGaugeLinks::U->N[2];z++){
		  for(INT y=0;y<=SUNcGaugeLinks::U->N[1]-1;y++){
			 for(INT x=0;x<=SUNcGaugeLinks::U->N[0]-1;x++){
				for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    
					JmOutStream  << x << " " << y << " " << z<<" "<< a;
					JmOutStream << " "<< SUNcChargeDensity::JpStat->Get(x,y,z,a)[0] << " "<< SUNcChargeDensity::Jp->Get(x,y,z,a)[0]<< " "<<SUNcChargeDensity::JmStat->Get(x,y,z,a)[0]<<" "<<SUNcChargeDensity::Jm->Get(x,y,z,a)[0] << " "<<std::endl;
				 }
			  }
		   }
	   }
	   
	   // CLOSE OUTPUT STREAM //
	   JmOutStream.close();
	}
    
    
    ////////////////////////////
    //     ENERGY PROFILE     //
    ////////////////////////////
	
	DOUBLE cB[Lattice::Dimension];
	DOUBLE cE[Lattice::Dimension];
	
	void SetConstants(){
		
		for(int mu=0;mu<Lattice::Dimension;mu++){
			
			cB[mu]=SQR(Lattice::a[mu]*SQR(Lattice::aScale)/Lattice::aCube);
			cE[mu]=SQR(Lattice::a[mu]*SQR(Lattice::aScale)/Lattice::aCube);
		}
	}

    void SaveFields(std::string EBfname,GaugeLinks *U){
		
	   std::ofstream EBOutStream;
	   
	   // OUTPUT FILES //
	   std::string EBOutputFile=StringManipulation::StringCast(IO::OutputDirectory,EBfname,"ID",RandomNumberGenerator::MySEED,".txt");
	   
	   // OPEN FILES //
	   EBOutStream.open(EBOutputFile.c_str());
	   
	   // SET PRECISION //
	   EBOutStream.precision(OUTPUT_PRECISION);
		
	   SET_ELEMENTARY_PLAQUETTE_BUFFERS();

	   // CREATE OUTPUT FOR DIAGONAL COMPONENT OF STRESS ENERGY TENSOR AND SQUARED FIELDS //
	   for(INT z=0;z<=SUNcGaugeLinks::U->N[2]-1;z++){
		   
		  DOUBLE EE0Sqr=0.0;   DOUBLE EE1Sqr=0.0;    DOUBLE EE2Sqr=0.0;
		  DOUBLE BB0Sqr=0.0;   DOUBLE BB1Sqr=0.0;    DOUBLE BB2Sqr=0.0;

		  for(INT y=0;y<=SUNcGaugeLinks::U->N[1]-1;y++){
			 for(INT x=0;x<=SUNcGaugeLinks::U->N[0]-1;x++){
				
					COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z,U);
					
					BB0Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uyz)/(Lattice::N[0]*Lattice::N[1]);
					BB1Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uzx)/(Lattice::N[0]*Lattice::N[1]);
					BB2Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uxy)/(Lattice::N[0]*Lattice::N[1]);
					
					for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    EE0Sqr+=SQR(SUNcElectricFields::E->Get(x,y,z,0,a)[0])/(Lattice::N[0]*Lattice::N[1]);
				    EE1Sqr+=SQR(SUNcElectricFields::E->Get(x,y,z,1,a)[0])/(Lattice::N[0]*Lattice::N[1]);
				    EE2Sqr+=SQR(SUNcElectricFields::E->Get(x,y,z,2,a)[0])/(Lattice::N[0]*Lattice::N[1]);

				}
			 }
		  }
		  
		
		  //OUTPUT TRANSVERSE AND LONGITUDINAL ENERGY DENSITY //
		   EBOutStream << z << "  " << DOUBLE(0.5)*(BB0Sqr+BB1Sqr)<< " "<<DOUBLE(0.5)*BB2Sqr<<"  "<<DOUBLE(0.5)*(EE0Sqr+EE1Sqr)<<"  "<<DOUBLE(0.5)*EE2Sqr;
		   EBOutStream << std::endl;
	   }
			
	   EBOutStream.close();
    }

	void SaveTransverselyIntegratedCurrents(std::string Jmfname){
	
	   std::ofstream JmOutStream;
	
	   // OUTPUT FILES //
	   std::string JmOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Jmfname,"ID",RandomNumberGenerator::MySEED,".txt");
	
	   // OPEN FILES //
	   JmOutStream.open(JmOutputFile.c_str());
	
	   // SET PRECISION //
	   JmOutStream.precision(OUTPUT_PRECISION);
	
	   // CREATE OUTPUT //
	   for(INT z=-1;z<=SUNcGaugeLinks::U->N[2];z++){
	
		  DOUBLE TrJmJm=0.0;  DOUBLE TrJpJp=0.0;
		  DOUBLE TrJmStJmSt=0.0;  DOUBLE TrJpStJpSt=0.0;
	
		  DOUBLE JmJm=0.0;   DOUBLE JpJp=0.0;
		  DOUBLE JmStJmSt=0.0;   DOUBLE JpStJpSt=0.0;
	
	
		  for(INT y=0;y<=SUNcGaugeLinks::U->N[1]-1;y++){
			 for(INT x=0;x<=SUNcGaugeLinks::U->N[0]-1;x++){
				for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
	
					JmJm+=SQR(SUNcChargeDensity::Jm->Get(x,y,z,a)[0]);
					JmStJmSt+=SQR(SUNcChargeDensity::JmStat->Get(x,y,z,a)[0]);
					JpJp+=SQR(SUNcChargeDensity::Jp->Get(x,y,z,a)[0]);
					JpStJpSt+=SQR(SUNcChargeDensity::JpStat->Get(x,y,z,a)[0]);
				}
			 }
		  }
	
		  TrJmJm=DOUBLE(0.5)*(JmJm)/(Lattice::N[0]*Lattice::N[1]);
		  TrJpJp=DOUBLE(0.5)*(JpJp)/(Lattice::N[0]*Lattice::N[1]);
		  TrJmStJmSt=DOUBLE(0.5)*(JmStJmSt)/(Lattice::N[0]*Lattice::N[1]);
		  TrJpStJpSt=DOUBLE(0.5)*(JpStJpSt)/(Lattice::N[0]*Lattice::N[1]);
	
		  // OUTPUT CHARGE DENSITY Jm //
		  JmOutStream << z <<"  " << TrJmJm << "  " << TrJmStJmSt << "  " << TrJpJp << "  " << TrJpStJpSt <<std::endl;
	   }
	
	   // CLOSE OUTPUT STREAM //
	   JmOutStream.close();
	}
}

#endif



