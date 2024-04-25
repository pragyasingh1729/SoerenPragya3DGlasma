//INCLUDE BASIC DEFINITIONS AND COMPILER FLAGS
#include "Definitions.cpp"

//SET GAUGE GROUP
#define SU_Nc_FLAG SU2_FLAG

// INCLUDE TIMING //
#include "MISC/Timing.cpp"

// INCLUDE SIMULATION PROCEDURE //
// #include "Simulation.cpp"

//INCLUDE SIMULATION PROCEDURE FOR PHYSICAL CHARGES //
// #include "SimulationmitPC.cpp"

//%%%%%%%%%%%%%%%%%%%%%%%%//
//TO BE CHANGED IN FUTURE//
//%%%%%%%%%%%%%%%%%%%%%%%//

// INCLUDE SIMULATION PROCEDURE FOR LEFT NUCLEUS //
// #include "LeftNucleusSimulation.cpp"

// INCLUDE SIMULATION PROCEDURE FOR RIGHT NUCLEUS //
#include "RightNucleusSimulation.cpp"

// GLOBAL RANDOM GENERATOR SEED //
INT GLOBAL_RNG_SEED;

// INCLUDE COMMANDLINE PROCESSING //
#include "IO/COMMANDLINE/cfile.c"
#include "CommandlineArguments.cpp"


//MAIN
int main(int argc,char **argv){
    
    ////////////////////////////////////
    // INITIALIZE OPEN MP ENVIRONMENT //
    ////////////////////////////////////
    
    //CHECK THREAD COUNT
    
    std::cerr << "#NUMBER OF THREADS " << omp_get_max_threads() << std::endl;
    
    //INITIALIZE THREADED FFTW
    int FFTW3_THREAD_STATUS=fftw_init_threads();
    std::cerr << "#FFTW THREAD STATUS " << FFTW3_THREAD_STATUS << std::endl;
    
    if(FFTW3_THREAD_STATUS==1){
        fftw_plan_with_nthreads(omp_get_max_threads());
    }
    
    // COMMANDLINE OUTPUT OF GAUGE GROUP //
    std::cerr << "#GAUGE GROUP IS SU(" << Nc << ")" << std::endl;
    
    // COMMANDLINE OUTPUT OF WORKING PRECISION //
    std::cerr << "#PRECISION IS " << MAX_DIGITS_PRECISION << " DIGITS" << std::endl;

    
    //////////////////////////////////
    //PROCESS COMMANDLINE ARGUMENTS //
    //////////////////////////////////
    
    ProcessCommandlineArguments(argc,argv);
	
    //////////////
    // SIMULATE //
    //////////////
    
    Timing::Reset();
	
	//INITIALIZE SIMULATION
	Simulation::Init();
	
	if(Parameters::Load==0){
	
		// COMMADNLINE NOTIFICATION //
		std::cerr << "#STARTING ID=" << GLOBAL_RNG_SEED << std::endl;
		
		// PERFORM CLASSICAL STATISTICAL SIMULATION //
		Simulation::Run(GLOBAL_RNG_SEED);
	}
	
	else{
		
		Simulation::Run(Dynamics::InputOne,Dynamics::InputTwo,Dynamics::InputThree,Dynamics::LoadTime);
	}
    // COMMADNLINE NOTIFICATION //
    std::cerr << "#COMPLETED ID=" << GLOBAL_RNG_SEED << std::endl;
	
    //EXIT
    exit(0);
}








