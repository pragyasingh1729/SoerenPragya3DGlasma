namespace InitialConditions{


    namespace SUNcChargeDensities{
        
        void SetRandom(DOUBLE JpAmpl,DOUBLE JmAmpl){
            
            std::cerr << "#SETTING RANDOM CURRENTS JP,JM" << std::endl;
            
            DOUBLE JpStatNetCharge[SUNcAlgebra::NumberOfGenerators];
            DOUBLE JmStatNetCharge[SUNcAlgebra::NumberOfGenerators];
            DOUBLE JpNetCharge[SUNcAlgebra::NumberOfGenerators];
            DOUBLE JmNetCharge[SUNcAlgebra::NumberOfGenerators];
            
            for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                JpStatNetCharge[a]=0.0;
                JmStatNetCharge[a]=0.0;
                JpNetCharge[a]=0.0;
                JmNetCharge[a]=0.0;
               
            }
            
            for(INT z=0;z<=Lattice::N[2]-1;z++){
                for(INT y=0;y<=Lattice::N[1]-1;y++){
                    for(INT x=0;x<=Lattice::N[0]-1;x++){
                        
                        
                        for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            
                            DOUBLE JpStat=JpAmpl*RandomNumberGenerator::Gauss();
                            DOUBLE JmStat=JmAmpl*RandomNumberGenerator::Gauss();

                            SUNcChargeDensity::Jp->Get(x,y,z,a)[0]=0;
                            SUNcChargeDensity::Jm->Get(x,y,z,a)[0]=0;
                            
                            SUNcChargeDensity::JpStat->Get(x,y,z,a)[0]=JpStat;
                            SUNcChargeDensity::JmStat->Get(x,y,z,a)[0]=JmStat;
                            
                            JpStatNetCharge[a]+=JpStat; JmStatNetCharge[a]+=JmStat;
                            //JpNetCharge[a]+=Jp;   JmNetCharge[a]+=Jm;
                        }
                    }
                }
            }
            
            for(INT z=0;z<=Lattice::N[2]-1;z++){
                for(INT y=0;y<=Lattice::N[1]-1;y++){
                    for(INT x=0;x<=Lattice::N[0]-1;x++){
                        
                        for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){                            
                           
                            SUNcChargeDensity::JpStat->Get(x,y,z,a)[0]-=JpStatNetCharge[a]/(Lattice::N[2]*Lattice::N[1]*Lattice::N[0]);
                            SUNcChargeDensity::JmStat->Get(x,y,z,a)[0]-=JmStatNetCharge[a]/(Lattice::N[2]*Lattice::N[1]*Lattice::N[0]);
                            
                        }
                        
                    }
                }
            }
            
            
            // CHANGE ELECTRIC FIELDS TO SATISFY GAUSS LAW //
            GaussLawRestoration::Restore();
            
            SU_Nc_ALGEBRA_FORMAT JpUp[SUNcAlgebra::NumberOfGenerators];
            SU_Nc_ALGEBRA_FORMAT JmUp[SUNcAlgebra::NumberOfGenerators];
            
           for(INT z=0;z<=Lattice::N[2]-1;z++){
                for(INT y=0;y<=Lattice::N[1]-1;y++){
                    for(INT x=0;x<=Lattice::N[0]-1;x++){
                        
                        
                        for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            
                            SUNcAlgebra::Operations::AdjointMultiplication(SUNcGaugeLinks::U->Get(x,y,z,2),SUNcChargeDensity::JpStat->Get(x,y,z+1,0),JpUp);
                            SUNcAlgebra::Operations::AdjointMultiplication(SUNcGaugeLinks::U->Get(x,y,z,2),SUNcChargeDensity::JmStat->Get(x,y,z+1,0),JmUp);
    
                            SUNcChargeDensity::Jp->Get(x,y,z,a)[0]=DOUBLE(0.5)*(SUNcChargeDensity::JpStat->GetValue(x,y,z,a)+JpUp[a]);
                            SUNcChargeDensity::Jm->Get(x,y,z,a)[0]=DOUBLE(0.5)*(SUNcChargeDensity::JmStat->GetValue(x,y,z,a)+JmUp[a]);
                            
                            
                        }
                    }
                }
            }
       
           // PERFORM INITIAL UPDATE STEP OF Jp AND Jm //
           Dynamics::SUNc::UpdateJ(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,SUNcGaugeLinks::U,SUNcChargeDensity::Jp,SUNcChargeDensity::Jm,SUNcChargeDensity::JpStat,SUNcChargeDensity::JmStat,DOUBLE(0.5)*Dynamics::dTStep);
		  
            
        }
    }
}
