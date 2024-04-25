namespace IO{
    
    void LoadGaugeLinks(std::string Ufname,GaugeLinks *U){
        
		// COMMANDLINE OUTPUT //
		std::cerr << "#LOADING GAUGE LINKS " << Ufname << std::endl;
		
		//LOAD FROM FILE
		std::ifstream UIn;
		
		UIn.open(Ufname.c_str());
		
		std::string ULink;
		
		//GET POSITION VALUES
		INT xIN,yIN,zIN,muIN;
		
		//NUMBER OF LINES READ FROM FILE
		INT InputCount=0;
		
		// MONITOR PRECISION //
		DOUBLE MaxUnitarityViolation=DOUBLE(0.0);
		
		//GET GAUGE LINK DATA FROM INPUT FILES
		while(UIn.good()){
			
			//READ LINES
			getline(UIn,ULink);
			
			//PROCESS FILE LINE BY LINE
			if(!(ULink.empty())){
				
				//STRING TOKEN
				std::stringstream ULinkValues(ULink);
				
				//GET POSITIONS IN FILE
				ULinkValues >> xIN; ULinkValues >> yIN; ULinkValues >> zIN;  ULinkValues >> muIN;

				//CHECK POSITIONS AND SET VALUES TO GAUGE LINK ARRAY
				if(InputCount==GaugeLinkIndex(xIN,yIN,zIN,muIN)){
					
					std::stringstream UString;    std::string strBuff;
					
					while(ULinkValues >> strBuff){UString << strBuff << " ";}
					
					SUNcGroup::IO::StringToMatrix(UString.str(),U->Get(xIN,yIN,zIN,muIN));
											
				}
				
				//BREAK IF POSITIONS DO NOT MATCH
				else{
					std::cerr << "#GAUGE LINK INPUT CAUSED FATAL ERROR -- TRANSVERSE POSITIONS DO NOT MATCH" << std::endl;
					std::cerr << InputCount << " " << xIN << " " << yIN << " " << zIN << " " << muIN << std::endl;
					exit(0);
				}
				
				//INCREASE POSITION COUNT
				InputCount+=SUNcGroup::MatrixSize;
			}
		}
				std::cout<<InputCount<<"  "<<SUNcGroup::MatrixSize*Lattice::Dimension*Lattice::N[0]*Lattice::N[1]*Lattice::N[2]<<std::endl;
		
				if(InputCount!=SUNcGroup::MatrixSize*Lattice::Dimension*Lattice::N[0]*Lattice::N[1]*(Lattice::N[2]+2)){
					std::cerr << "#GAUGE LINK INPUT CAUSED FATAL ERROR -- GAUGE LINKS NOT LOADED CORRECTLY" << std::endl;
					exit(0);
				}
		
				std::cerr << "#MAX UNITARITY VIOLATION IS " << MaxUnitarityViolation << std::endl;
    }
    
    
    void LoadElectricFields(std::string Efname,ElectricFields *E){
        
            // COMMANDLINE OUTPUT //
            std::cerr << "#LOADING ELECRIC FIELDS " << Efname << std::endl;
            
            //LOAD FROM FILE
            std::ifstream EIn;
            
            EIn.open(Efname.c_str());
            
            std::string EAlgeb;
            
            //GET POSITION VALUES
            INT xIN,yIN,zIN,muIN;
            
            //NUMBER OF LINES READ FROM FILE
            INT InputCount=0;
            
            //GET GAUGE LINK DATA FROM INPUT FILES
            while(EIn.good()){
                
                //READ LINES
                getline(EIn,EAlgeb);
                
                //PROCESS FILE LINE BY LINE
                if(!(EAlgeb.empty())){
                    
                    //STRING TOKEN
                    std::stringstream EAlgebValues(EAlgeb);
                    
                    //GET POSITIONS IN FILE
                    EAlgebValues >> xIN; EAlgebValues >> yIN; EAlgebValues >> zIN; EAlgebValues >> muIN;
					
                    //CHECK POSITIONS AND SET VALUES TO GAUGE LINK ARRAY
                    if(InputCount==VectorFieldIndex(xIN,yIN,zIN,muIN,0)){
                        
                        DOUBLE CurrentValue;
                        
                        for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            
                            if(EAlgebValues.good()){
                                EAlgebValues >> CurrentValue;
                            }
                            else{
                                std::cerr << "#COULD NOT READ COLOR VECTOR CORRECTLY" << std::endl;
                                exit(0);
                            }
                            
                            E->Get(xIN,yIN,zIN,muIN,a)[0]=CurrentValue;
                            
                        }
                        
                    }
                    
                    //BREAK IF POSITIONS DO NOT MATCH
                    else{
                        std::cerr << "#ELECTRIC FIELD INPUT CAUSED FATAL ERROR -- TRANSVERSE POSITIONS DO NOT MATCH" << std::endl;
                        exit(0);
                    }
                    
                    //INCREASE POSITION COUNT
                    InputCount+=SUNcAlgebra::NumberOfGenerators;
                }
            }
            
					std::cout<<InputCount<<"  "<<SUNcAlgebra::NumberOfGenerators*Lattice::Dimension*Lattice::N[0]*Lattice::N[1]*Lattice::N[2]<<std::endl;

				if(InputCount!=SUNcAlgebra::NumberOfGenerators*Lattice::Dimension*Lattice::N[0]*Lattice::N[1]*(Lattice::N[2]+2)){
					std::cerr << "#ELECTRIC FIELD INPUT CAUSED FATAL ERROR -- ELECTRIC FIELDS NOT LOADED CORRECTLY" << std::endl;
					exit(0);
				}
    }
    
	 void LoadCurrents(std::string Jfname,ChargeDensity *JpStat, ChargeDensity *Jp, ChargeDensity *JmStat, ChargeDensity *Jm){

		// COMMANDLINE OUTPUT //
		std::cerr << "#LOADING CURRENTS " << Jfname << std::endl;

		//LOAD FROM FILE
		std::ifstream JIn;

		JIn.open(Jfname.c_str());

		std::string JAlgeb;

		//GET POSITION VALUES
		INT xIN,yIN,zIN,aIn;

		//GET VALUES
		DOUBLE JpStatIN, JpIN,JmIN,JmStatIN;

		//GET GAUGE LINK DATA FROM INPUT FILES
		while(JIn.good()){

			//READ LINES
			getline(JIn,JAlgeb);

			//PROCESS FILE LINE BY LINE
			if(!(JAlgeb.empty())){

				//STRING TOKEN
				std::stringstream JAlgebValues(JAlgeb);

				//GET POSITIONS IN FILE
				JAlgebValues >> xIN; JAlgebValues >> yIN; JAlgebValues >> zIN;  JAlgebValues>>aIn; JAlgebValues >> JpStatIN; JAlgebValues >> JpIN;  JAlgebValues >> JmStatIN;  JAlgebValues >> JmIN;

				
				for(aIn==0;aIn<3;aIn++){
					JpStat->Get(xIN,yIN,zIN,aIn)[0]=JpStatIN;
					Jp->Get(xIN,yIN,zIN,aIn)[0]=JpIN;
					JmStat->Get(xIN,yIN,zIN,aIn)[0]=JmStatIN;
					Jm->Get(xIN,yIN,zIN,aIn)[0]=JmIN;

				}
			}
		}
	}

   
	void LoadConfiguration(std::string Ufname,std::string Efname,std::string Jfname,GaugeLinks *U,ElectricFields *E,ChargeDensity *JpStat,ChargeDensity *Jp,ChargeDensity *JmStat,ChargeDensity *Jm){
        
        LoadGaugeLinks(Ufname,U);
        LoadElectricFields(Efname,E);
		LoadCurrents(Jfname,JpStat,Jp,JmStat,Jm);
    }
    
    void LoadConfiguration(std::string Ufname,std::string Efname,std::string Jfname){
        
		LoadConfiguration(Ufname,Efname,Jfname,SUNcGaugeLinks::U,SUNcElectricFields::E,SUNcChargeDensity::JpStat,SUNcChargeDensity::Jp,SUNcChargeDensity::JmStat,SUNcChargeDensity::Jm);
        
    }
    
}







