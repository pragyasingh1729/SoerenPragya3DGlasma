#ifndef _COMMANDLINE_ARGUMENTS_CYM_SIMULATION_
#define _COMMANDLINE_ARGUMENTS_CYM_SIMULATION_

void ProcessCommandlineArguments(int argc,char **argv){
    
    Konfig arguments(argc,argv);
    
    ////////////////////////////////
    // SET INPUT/OUTPUT DIRECTORY //
    ////////////////////////////////
    
    char OutDir[512]="OUTPUT";
    
    arguments.Getval("o",OutDir);
    
    IO::SetOutputDirectory(OutDir);
    
    ////////////////////////////
    // DETERMINE LATTICE SIZE //
    ////////////////////////////
    
    INT NsSites=-1;
    
    arguments.Getval("N",NsSites);
    
    INT NxSites=NsSites;
    INT NySites=NsSites;
    INT NzSites=NsSites;

    arguments.Getval("Nx",NxSites);
    arguments.Getval("Ny",NySites);
    arguments.Getval("Nz",NzSites);

    if((NxSites>0)&&(NySites>0)&&(NzSites>0)){

        Lattice::N[0]=NxSites;
        Lattice::N[1]=NySites;
        Lattice::N[2]=NzSites;

        Lattice::Volume=NxSites*NySites*NzSites;

        std::cerr << "## LATTICE SIZE IS " << Lattice::N[0] << "x" << Lattice::N[1] << "x" << Lattice::N[2] <<  std::endl;
    }
    else{
        std::cerr << "## NUMBER OF LATTICE SITES NOT SPECIFIED -- USING " << Lattice::N[0] << "x" << Lattice::N[1] << "x" << Lattice::N[2]  << std::endl;
    }
    
    
    DOUBLE aSpacing=-1.0;
    
    arguments.Getval("a",aSpacing);
    
    DOUBLE aXSpacing=aSpacing;
    DOUBLE aYSpacing=aSpacing;
    DOUBLE aZSpacing=aSpacing;
    
    arguments.Getval("ax",aXSpacing);
    arguments.Getval("ay",aYSpacing);
    arguments.Getval("az",aZSpacing);
    
    if((aXSpacing>0.0)&&(aYSpacing>0.0)&&(aZSpacing>0.0)){
        
        Lattice::a[0]=aXSpacing;
        Lattice::a[1]=aYSpacing;
        Lattice::a[2]=aZSpacing;
        
        Lattice::aCube=aXSpacing*aYSpacing*aZSpacing;
        
        std::cerr << "## LATTICE SPACINGS ARE " << Lattice::a[0] << " " << Lattice::a[1] << " " << Lattice::a[2] <<  std::endl;
    }
    else{
        std::cerr << "## LATTICE SPACING NOT SPECIFIED -- USING " << Lattice::a[0] << " " << Lattice::a[1] << " " << Lattice::a[2] <<  std::endl;
    }

    
    arguments.Getval("dT",Dynamics::dTStep);
    
    std::cerr << "## DYNAMICAL TIME STEP IS " << Dynamics::dTStep <<  std::endl;
    
    arguments.Getval("RlOverAzL",Parameters::SigmaLOverAzL);
	
	arguments.Getval("RlOverAzR",Parameters::SigmaLOverAzR);
	
	arguments.Getval("CoMEnergy",Parameters::CoMEnergy);
	
    arguments.Getval("QsA",Parameters::Qs0A);
		
    arguments.Getval("MEff",Parameters::MEff);
	
    arguments.Getval("Lambda",Parameters::Lambda);
        
    arguments.Getval("SEED",GLOBAL_RNG_SEED);
	
	arguments.Getval("i1",Dynamics::InputOne);

	arguments.Getval("i2",Dynamics::InputTwo);

	arguments.Getval("i3",Dynamics::InputThree);

	arguments.Getval("Time",Dynamics::LoadTime);
	
	arguments.Getval("Loaded",Parameters::Load);
	
	arguments.Getval("NucleusL",Parameters::LeftNucleus);
	
	arguments.Getval("NucleusR",Parameters::RightNucleus);
	
}
        
#endif
