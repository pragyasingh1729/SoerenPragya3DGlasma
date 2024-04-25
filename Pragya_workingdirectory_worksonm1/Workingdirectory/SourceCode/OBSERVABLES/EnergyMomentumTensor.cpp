#ifndef __ENERGYDENSITYMAP__CPP__
#define __ENERGYDENSITYMAP__CPP__

namespace Observables {
        
    namespace EnergyMomentumTensor{
		
		//COMPUTES COMPONENT OF ENERGY MOMENTUM TENSOR BY AVERAGING OVER LONGITUDINAL PLANE FOR EVENT MONITOR //

        void CreateMapBfCollision(std::string fname,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
            
            //CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
            DOUBLE cB[Lattice::Dimension];  DOUBLE cE[Lattice::Dimension];
            
            //SET CONSTANTS
            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                cB[mu]=SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                cE[mu]=SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
            }
			
            //BUFFERS FOR AVERAGE FIELD STRENGTH
            SET_AVG_FIELD_STRENGTH_BUFFERS();
            
            //CREATE OUTPUT STREAM
            std::ofstream OutStream;
            
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            OutStream.open(OutputFile.c_str());
			
			OutStream.precision(OUTPUT_PRECISION);
			
			for(INT y=yLow;y<=yHigh;y++){
				for(INT x=xLow;x<=xHigh;x++){
					
					DOUBLE T00Left=0.0;		DOUBLE T00Right=0.0;

					for(INT z=zLow; z<=zHigh; z++){
						
						DOUBLE EE0Sqr=0.0;   DOUBLE EE1Sqr=0.0;    DOUBLE EE2Sqr=0.0;
						DOUBLE BB0Sqr=0.0;   DOUBLE BB1Sqr=0.0;    DOUBLE BB2Sqr=0.0;
						
						//CALCULATES E'S FOR SLICE t+dt/2 WHILE B'S FOR t+dt
						COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z,U);
						
						//ADD SQUARED FIELD STRENGTH TO LOCAL OBSERVABLES
						BB0Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uyz);
						BB1Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uzx);
						BB2Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uxy);

						for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
		
							EE0Sqr+=SQR(E->Get(x,y,z,0,a)[0]);
							EE1Sqr+=SQR(E->Get(x,y,z,1,a)[0]);
							EE2Sqr+=SQR(E->Get(x,y,z,2,a)[0]);
							
						}
					
						if(z<=(zHigh-1)/2){
							T00Left+=DOUBLE(0.5)*(cE[0]*EE0Sqr+cE[1]*EE1Sqr+cE[2]*EE2Sqr+cB[0]*BB0Sqr+cB[1]*BB1Sqr+cB[2]*BB2Sqr)/Lattice::N[2];

						}
						else{
							T00Right+=DOUBLE(0.5)*(cE[0]*EE0Sqr+cE[1]*EE1Sqr+cE[2]*EE2Sqr+cB[0]*BB0Sqr+cB[1]*BB1Sqr+cB[2]*BB2Sqr)/Lattice::N[2];
						}
					}
					
					OutStream << x << "  " << y << "  " << T00Left << "  "<< T00Right << std::endl;

				}
			}
			
		 //CLOSE OUTPUT STREAM
		 OutStream.close();
		}
		
        //CREATE ENERGY DENSITY MAP
        void CreateMapBfCollision(std::string fname,GaugeLinks *U,ElectricFields *E){
            
            CreateMapBfCollision(fname,0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
            
        }
        
        
        void CreateMapBfCollision(std::string fname){
            
            CreateMapBfCollision(fname,SUNcGaugeLinks::U,SUNcElectricFields::E);
            
        }
		
		void CreateMap(std::string fname,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,GaugeLinks *UOld,ElectricFields *E){
			
			//CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
			DOUBLE qB[Lattice::Dimension];  DOUBLE qE[Lattice::Dimension];
			
			//SET CONSTANTS
			for(int mu=0;mu<Lattice::Dimension;mu++){
				
				qB[mu]=(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
				qE[mu]=(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
			}
			
			
			//CREATE OUTPUT STREAM
			std::ofstream OutStream;

			OutStream.precision(OUTPUT_PRECISION);

			std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
			
			OutStream.open(OutputFile.c_str());
						
			//UPDATE ALL GAUGE LINKS
			for(INT z=zLow;z<=zHigh;z++){
				
				//DIAGONAL COMPONENTS OF THE ENERGY MOMENTUM TENSOR
				DOUBLE T00=0.0,TXX=0.0,TYY=0.0,TZZ=0.0;
				DOUBLE T0X=0.0,T0Y=0.0,T0Z=0.0,TXY=0.0,TYZ=0.0,TZX=0.0;
				
			#pragma omp parallel
				{
				//BUFFERS FOR AVERAGE FIELD STRENGTH
				SET_AVG_FIELD_STRENGTH_BUFFERS();

			#pragma omp for reduction( + : T00,TXX,TYY,TZZ,T0X,T0Y,T0Z,TXY,TYZ,TZX)

				for(INT y=yLow;y<=yHigh;y++){
					for(INT x=xLow;x<=xHigh;x++){
						
						//COMPUTE AVERAGE FIELD STRENGTH
						SUNcGaugeLinks::UOld->SynchronizeGhostCells();
						SUNcGaugeLinks::U->SynchronizeGhostCells();

						//CALCULATES E'S AND B'S FOR SLICE t+dt/2
						COMPUTE_AVG_FIELD_STRENGTH(x,y,z,UOld);
						
     					//COMPUTE DIAGONAL COMPONENTS OF ENERGY MOMENTUM TENSOR
						T00+=(DOUBLE(0.5)*(SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc+SQR(qB[2])*B2SqrLoc))/(Lattice::N[0]*Lattice::N[1]);
						TXX+=(DOUBLE(0.5)*((SQR(qE[1])*E1SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[1])*B1SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qE[0])*E0SqrLoc+SQR(qB[0])*B0SqrLoc)))/(Lattice::N[0]*Lattice::N[1]);
						TYY+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qE[1])*E1SqrLoc+SQR(qB[1])*B1SqrLoc)))/(Lattice::N[0]*Lattice::N[1]);
						TZZ+=(DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc)-(SQR(qE[2])*E2SqrLoc+SQR(qB[2])*B2SqrLoc)))/(Lattice::N[0]*Lattice::N[1]);
						
						//COMPUTE POYNTING VECTOR
						T0X+=(qE[1]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B2Loc)-qE[2]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B1Loc))/(Lattice::N[0]*Lattice::N[1]);
						T0Y+=(qE[2]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B0Loc)-qE[0]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B2Loc))/(Lattice::N[0]*Lattice::N[1]);
						T0Z+=(qE[0]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B1Loc)-qE[1]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B0Loc))/(Lattice::N[0]*Lattice::N[1]);

						//COMPUTE SPATIAL STRESS
						TXY+=(qE[0]*qE[1]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,E1Loc)+qB[0]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(B0Loc,B1Loc))/(Lattice::N[0]*Lattice::N[1]);
						TYZ+=(qE[1]*qE[2]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,E2Loc)+qB[1]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(B1Loc,B2Loc))/(Lattice::N[0]*Lattice::N[1]);
						TZX+=(qE[2]*qE[0]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,E0Loc)+qB[2]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(B2Loc,B0Loc))/(Lattice::N[0]*Lattice::N[1]);

				
					}
				}
			}
				//WRITE OUTPUT TO FILE
				OutStream << z   <<"  "<< T00<<"  "<< TXX << "  " << TYY<<"  "<<TZZ;
				OutStream << "  "<< T0X << "  "<< T0Y <<"  "<< T0Z;
				OutStream << "  "<< TXY << "  " <<TYZ << "  " << TZX;
			
				OutStream<<std::endl;
			
			}
			
			//CLOSE OUTPUT STREAM
			OutStream.close();
			
		}
		
		
		//CREATE ENERGY DENSITY MAP
		void CreateMap(std::string fname,GaugeLinks *U,GaugeLinks *UOld,ElectricFields *E){
			
			CreateMap(fname,0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,UOld,E);
			
		}
		
		void CreateMap(std::string fname){
			
			CreateMap(fname,SUNcGaugeLinks::U,SUNcGaugeLinks::UOld,SUNcElectricFields::E);
			
		}
		
    }
}
#endif
