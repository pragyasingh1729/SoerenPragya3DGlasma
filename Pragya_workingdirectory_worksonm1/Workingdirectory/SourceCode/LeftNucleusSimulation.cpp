/////////////////////////////////////////
// SU(Nc) GROUP AND ALGEBRA OPERATIONS //
/////////////////////////////////////////

//INCLUDE SU(Nc) ALGEBRA AND GROUP DEFINITION
#include "SUNc/Definitions.cpp"

/////////////////////////////
// RANDOM NUMBER GENERATOR //
/////////////////////////////

//INCLUDE RANDOM NUMBER GENERATOR
#include "MISC/RNG/GSLRandomNumberGenerator.cpp"

/////////////////////////////////////////////
// DEFINITION OF LATTICE GRID AND INDEXING //
/////////////////////////////////////////////

#include "LATTICE/3DGrid.cpp"
#include "LATTICE/Parameters.cpp"
#include "LATTICE/WilsonLines.cpp"
#include "LATTICE/GaugeLinks.cpp"
#include "LATTICE/ElectricFields.cpp"
#include "LATTICE/GaugeTransformations.cpp"
#include "LATTICE/ChargeDensity.cpp"
#include "LATTICE/SUNcVariables.cpp"
#include "LATTICE/Indexing.cpp"
#include "LATTICE/FourierSpace.cpp"


///////////////////////////
//GENERAL PURPOSE MACROS //
///////////////////////////

//INCLUDE MACROS TO COMPUTE PLAQUETTES
#include "CALC/Plaquettes.cpp"

//INCLUDE MACROS TO DETERMINE AVERAGE FIELD STRENGTH
#include "CALC/AvgFieldStrength.cpp"

//INCLUDE MACROS TO COMPUTE GAUSS LAW VIOLATION
#include "CALC/GaussViolation.cpp"

//INCLUDE MACROS TO PERFORM GAUGE TRANSFORMATIONS
#include "CALC/GaugeTransformation.cpp"

#include "SUNcIndexing.cpp"

///////////////////////////
//  OUTPUT HANDLING      //
///////////////////////////

#include "IO/StringManipulation.cpp"
#include "IO/OutputManagement.cpp"
#include "IO/SaveConfiguration.cpp"
#include "IO/LoadConfiguration.cpp"


///////////////////////
//REAL-TIME DYNAMICS //
///////////////////////

//INCLUDE BASIC DYNAMICS DEFINITIONS
#include "DYNAMICS/SUNcDynamics.cpp"

//////////////////////
//BASIC OBSERVABLES //
//////////////////////

//INCLUDE MEASUREMENT OF BULK OBSERVABLES
#include "OBSERVABLES/BulkObservables.cpp"

//INCLUDE MEASUREMENT OF GAUSS LAW VIOLATION
#include "OBSERVABLES/GaussLawViolation.cpp"

//INCLUDE MEASUREMENT OF UNITARITY VIOLATION
#include "OBSERVABLES/UnitarityViolation.cpp"

//INCLUDE MEASUREMENT OF ENERGY MOMENTUM TENSOR
#include "OBSERVABLES/EnergyMomentumTensor.cpp"

////////////////////////////
//GAUGE FIXING PROCEDURES //
////////////////////////////

//INCLUDE GAUGE TRANSFORMATION RULES AND BUFFERS
#include "GAUGETRANSFORMATION/GaugeTransformation.cpp"

//INCLUDE COULOMB GAUGE FIXING PROCEDURE
#include "GAUGETRANSFORMATION/CoulombGaugeFixing.cpp"

///////////////////////
//INITIAL CONDITIONS //
///////////////////////

#include "FFT/FFT1D.cpp"
#include "FFT/FFT2D.cpp"
#include "FFT/FFT3D.cpp"

#include "COLLIDINGOBJECTS/CollidingObjects.cpp"
#include "COLLIDINGOBJECTS/Comparision.cpp"

#include "INITIALCONDITIONS/SUNCCURRENTS/RevSingleNucleus.cpp"

///////////////////////////
// SIMULATION PROCEDURE  //
///////////////////////////

//SIMULATION PARAMETERS

namespace Simulation {
   
    INT MY_MPI_RNG_SEED;
    
    void Init(){
	   
	   /////////////////////
	   // INITIALIZATIONS //
	   /////////////////////
	   
	   // INITIALIZE 3D LATTICE SETUP //
	   Lattice::Init();
	   
	   // INITIALIZE GAUGE TRANSFORMATIONS //
	   GaugeTransformation::Init();
	   
	   // INITIALIZE GAUGE TRANSFORMED FIELDS //
	   GaugeFixedVariables::Init();
		
    }
    
    /////////////////////
    //CREATE INFO FILE //
    /////////////////////
    
    void CreateInfoFile(){
	   
	   // CREATE INFO FILE //
	   std::ofstream OutStream;
	   OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"InfoID",MY_MPI_RNG_SEED,".txt").c_str());
	   
	   // CREATE INFO FILE CONTAINING PARAMETERS //
	   OutStream << "#LATTICE DATA" << std::endl;
	   OutStream << "#Nx,Ny,Nz= " << SUNcGaugeLinks::U->N[0] << " " << SUNcGaugeLinks::U->N[1] << " " << SUNcGaugeLinks::U->N[2] << std::endl;
	   OutStream << "#ax,ay,az= " << SUNcGaugeLinks::U->a[0] << " " << SUNcGaugeLinks::U->a[1] << " " << SUNcGaugeLinks::U->a[2] << std::endl;
	   OutStream << "#dT= " << Dynamics::dTStep << std::endl;
	   OutStream << "    " << std::endl;
	   OutStream << "QsR= " << Parameters::Qs0R << std::endl;
	   OutStream << "QsAperp= " << Parameters::Qs0A << std::endl;
	   OutStream << "MEff aka Infrared Cutoff= " << double(Parameters::MEff/Parameters::Qs0A) << " "<<"(times QsAperp)"<<std::endl;
	   OutStream << "Lambda aka UV Cutoff= " << double(Parameters::Lambda/Parameters::Qs0A) << " "<<"(times QsAperp)"<<std::endl;
	   OutStream << "SigmaLOverAz L= " << Parameters::SigmaLOverAzL << std::endl;
	   OutStream << "SigmaLOverAz R= " << Parameters::SigmaLOverAzR << std::endl;
	   OutStream << std::endl;
	   
	   // CLOSE OUPUT STREAM //
	   OutStream.close();
	   
    }
    
	   //////////////////////////////////////////////
	   // CHECKS: WILSON LINES FOR LEFT AND RIGHT //
	   ////////////////////////////////////////////
    
    void WLinesLeftNucleus(){
	   
	   std::ofstream LSitOutStream;
	   LSitOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"LeftWilson",".txt").c_str());
	   
	   // SET PRECISION //
	   LSitOutStream.precision(OUTPUT_PRECISION);
	   
	   // CREATE OUTPUT //
	   for(INT y=0;y<=Lattice::N[1]-1;y++){
		  for(INT x=0;x<=Lattice::N[0]-1;x++){
			 
			 // OUTPUT WILSON LINES //
			 LSitOutStream << x << " " << y << " " << SUNcGroup::IO::MatrixToString(SUNcWilsonLines::VLeft->Get(x,y,0))<< std::endl;
		  }
	   }
	   
	   //CLOSE OUTPUT STREAM
	   LSitOutStream.close();
    }
    
    void WLinesRightNucleus(){
	   
	   std::ofstream RSitOutStream;
	   RSitOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"RightWilson",".txt").c_str());
	   
	   // SET PRECISION //
	   RSitOutStream.precision(OUTPUT_PRECISION);
	   
	   // CREATE OUTPUT //
	   for(INT y=0;y<=Lattice::N[1]-1;y++){
		  for(INT x=0;x<=Lattice::N[0]-1;x++){
			 
			 // OUTPUT WILSON LINES //
			 RSitOutStream << x << " " << y << " " << SUNcGroup::IO::MatrixToString(SUNcWilsonLines::VRight->Get(x,y,Lattice::N[2]-1))<< std::endl;
		  }
	   }
	   
	   //CLOSE OUTPUT STREAM
	   RSitOutStream.close();
    }
    
    
    ////////////////////////////////////
    // CLASSICAL YANG-MILLS EVOLUTION //
    ////////////////////////////////////
    
    void Evolve(DOUBLE MaxTime){

	   ///////////////////////
	   // OUTPUT MANAGEMENT //
	   ///////////////////////
	   
	   std::ofstream EnergyOutStream;
	   EnergyOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"EnergyID",MY_MPI_RNG_SEED,".txt").c_str());
	   EnergyOutStream.precision(OUTPUT_PRECISION);
	   
	   ///////////////////////////////
	   // INITIAL STATE OBSERVABLES //
	   ///////////////////////////////
	   
	   std::cerr << "#COMPUTING ENERGY DENSITY AT T=" << Dynamics::Time() << " ELAPSED TIME " << Timing::Get() << " s" << std::endl;
	
	   Timing::Reset();
	   
	   // COMPUTE BULK OBSERVABLES FOR GAUGE SECTOR //
	   Observables::Bulk::Update();
	   
	   // CREATE OUPUT //
	   EnergyOutStream << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX()<< " " << Observables::Bulk::TYY()<< " " << Observables::Bulk::TZZ()<< " " << Observables::Bulk::SOURCE() << std::endl;
	   
	   std::cerr << "#ENERGY DENSITY COMPUTED IN " << Timing::Get() << " s" << std::endl;
	   
	   std::cerr << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX()<< " " << Observables::Bulk::TYY()<< " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::SOURCE() << std::endl;
	   
	   Timing::Reset();
		
	   
	   ////////////////////
	   // TIME EVOLUTION //
	   ////////////////////
	   
	   // COMMANDLINE OUTPUT //
	   std::cerr << "#STARTING TIME EVOLUTION" << std::endl;
	   
	   // PERFORM EVOLUTION //
	   while(Dynamics::Time()<MaxTime){
		  
		  Dynamics::SUNc::Update(SUNcGaugeLinks::U,SUNcElectricFields::E,SUNcChargeDensity::Jp,SUNcChargeDensity::Jm,SUNcChargeDensity::JpStat,SUNcChargeDensity::JmStat);
		  
		  // EVOLUTION TIMING //
		  if(Dynamics::tSteps%10==0){
			 std::cerr << "#ELAPSED TIME " << Timing::Get() << " s FOR 10 EVOLUTION STEPS" << std::endl;
			 Timing::Reset();
		  }
		  
		  // COMPUTE BULK OBSERVABLES //
		  if(Dynamics::tSteps%20==0){
			 
			 std::cerr << "#COMPUTING ENERGY DENSITY AT T=" << Dynamics::Time() << std::endl;
			 Timing::Reset();
			 
			 // COMPUTE BULK OBSERVABLES FOR GAUGE SECTOR //
			 Observables::Bulk::Update();
			 
			 // CREATE OUPUT //
			 std::cerr << "#ENERGY DENSITY COMPUTED IN " << Timing::Get() << "s" << std::endl;
			 std::cerr << Dynamics::tSteps<< "  " <<Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::SOURCE() << std::endl;
			 
			 Timing::Reset();
		  }
		  
		  
		  if(Dynamics::tSteps%20==0){
			 EnergyOutStream <<Dynamics::tSteps<< "  " << Dynamics::Time() << "  " << Observables::Bulk::T00() << "  " << Observables::Bulk::SOURCE() << std::endl;
			 Observables::GaussLaw::SaveViolation(StringManipulation::StringCast("GausCheck"));
		  }
		  
		 //SAVE SU(Nc) GAUGE FIELD CONFIGURATIONS //
		  if(Dynamics::tSteps%100==0){
			 
			  std::cerr << "#SAVING OBSERVABLES AT T=" << Dynamics::tSteps<< "  " << Dynamics::Time() << std::endl;
			  Observables::EnergyMomentumTensor::CreateMap(StringManipulation::StringCast("EMTDynT",Dynamics::Time()));
		  }
		
			if(Dynamics::tSteps%3000==0){
			
				std::cerr << "#SAVING SU(Nc) GAUGE FIELD CONFIGURATIONS AT T=" << Dynamics::tSteps<< "  " << Dynamics::Time() << std::endl;
				
//				IO::SaveConfigurationJ(StringManipulation::StringCast("JDynT",Dynamics::Time()));
//				IO::SaveConfiguration(StringManipulation::StringCast("UDynT",Dynamics::Time()),StringManipulation::StringCast("EDynT",Dynamics::Time()));

			}
		  
		   // PRECISION CHECKS //
		  if(Dynamics::tSteps%1000==0){

			 std::cerr << "#CHECKING GAUSS LAW VIOLATION AND UNITARITY AT T=" << Dynamics::Time() << std::endl;
			 Timing::Reset();

			//CHECK GAUSS VIOLATION//
//			 Observables::GaussLaw::CheckViolation();
//			 Observables::GaussLaw::CheckViolation(StringManipulation::StringCast("GVDynT",Dynamics::Time()));

			 //CHECK UNITARITY VIOLATION//
			  
			 Observables::Unitarity::CheckViolation();

			 std::cerr << "#CHECKED IN " << Timing::Get() << " s" << std::endl;

			 Timing::Reset();
		  }
	   }
	   // CLOSE OUTPUT STREAMS //
	   EnergyOutStream.close();
	   
    }
    
    //////////////////////////
    //SIMULATION PROCDEDURE //
    //////////////////////////
    
    void Run(INT MPI_RNG_SEED){
	   
	   ///////////
	   // SETUP //
	   ///////////
	   
	   //SET SEED //
	   MY_MPI_RNG_SEED=MPI_RNG_SEED;
	   
	   //INITIALIZE RANDOM NUMBER GENERATOR //
	   RandomNumberGenerator::Init(MY_MPI_RNG_SEED);
	   
	   //INITIALIZE DYNAMICS //
	   Dynamics::Reset();
		
	   //SET LEFT SITTING NUCLEUS //
	   CollidingObjects::SetWilsonLines(SUNcChargeDensity::rhoOldL,SUNcChargeDensity::rhoL,SUNcWilsonLines::VLeft,'L',0.85*Lattice::N[2]/2-1, Parameters::Qs0A_L);

	   //SET RIGHT SITTING NUCLEUS //
	   CollidingObjects::SetWilsonLines(SUNcChargeDensity::rhoOldR,SUNcChargeDensity::rhoR,SUNcWilsonLines::VRight,'R',1.15*Lattice::N[2]/2-1, Parameters::Qs0A_L);
	
	   //////////////////////
	   // CREATE INFO FILE //
	   //////////////////////
	   CreateInfoFile();
	   
	   ////////////////////////////
	   // SET INITIAL CONDITIONS //
	   ////////////////////////////
	   
	   InitialConditions::SingleSheet::SingleSheet(SUNcChargeDensity::rhoL,0,Lattice::N[2]-1,SUNcChargeDensity::RhoLT1,SUNcWilsonLines::VLeft,'L',SUNcChargeDensity::JpStat,SUNcChargeDensity::Jp);
		
	   // SYNCHRONIZE FIELDS AGAIN //
		
	   SUNcChargeDensity::JpStat->SynchronizeGhostCells();
	   SUNcChargeDensity::Jp->SynchronizeGhostCells();
	   
	   SUNcChargeDensity::JmStat->SynchronizeGhostCells();
	   SUNcChargeDensity::Jm->SynchronizeGhostCells();
	   
	   SUNcGaugeLinks::U->SynchronizeGhostCells();
	   SUNcGaugeLinks::UOld->SynchronizeGhostCells();

	   SUNcElectricFields::E->SynchronizeGhostCells();
	   					   
	   Observables::GaussLaw::CheckViolation(StringManipulation::StringCast("GVDynT",Dynamics::Time()));
		
	   Observables::Unitarity::CheckViolation();
	 
		/////////////////////////////////////////////////////////
	   // SAVE INITIAL CONDITIONS //
	   /////////////////////////////////////////////////////////
	   
	  std::cerr << "#SAVING SU(Nc) GAUGE FIELD CONFIGURATIONS AT T=" << Dynamics::Time() << std::endl;
	  Observables::EnergyMomentumTensor::CreateMap(StringManipulation::StringCast("EMTDynT",Dynamics::Time()));

//	  IO::SaveConfigurationJ(StringManipulation::StringCast("JDynT",Dynamics::Time()));
//	  IO::SaveConfiguration(StringManipulation::StringCast("UDynT",Dynamics::Time()),StringManipulation::StringCast("EDynT",Dynamics::Time()));
				
	 
	   Observables::GaussLaw::SaveViolation(StringManipulation::StringCast("GausCheck"));
		
	   ////////////////////
	   // EVOLVE	     //
	   //////////////////
	   
	   Evolve(1100.0);
    }

	void Run(std::string f1name,std::string f2name,std::string f3name, DOUBLE t){
		
		IO::LoadConfiguration(f1name,f2name,f3name);
		
		SUNcChargeDensity::JpStat->SynchronizeGhostCells();
		SUNcChargeDensity::Jp->SynchronizeGhostCells();
		
		SUNcChargeDensity::JmStat->SynchronizeGhostCells();
		SUNcChargeDensity::Jm->SynchronizeGhostCells();
		
		SUNcGaugeLinks::U->SynchronizeGhostCells();
		SUNcElectricFields::E->SynchronizeGhostCells();
						
		Dynamics::LoadedTime(t);
		
		IO::SaveConfigurationJ(StringManipulation::StringCast("JDynT",Dynamics::Time()));
		IO::SaveConfiguration(StringManipulation::StringCast("UDynT",Dynamics::Time()),StringManipulation::StringCast("EDynT",Dynamics::Time()));
		
		Evolve(100.0);
	
	}
}





