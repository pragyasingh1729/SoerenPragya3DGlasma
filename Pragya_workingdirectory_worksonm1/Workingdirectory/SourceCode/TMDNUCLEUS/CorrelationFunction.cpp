#ifndef __CHARGE_DISTRIBUTIONS__
#define __CHARGE_DISTRIBUTIONS__

namespace ChargeDistributions{
	
	FFT3D *FourierVariables;
		
	INT SUNcAlgebraIndex3D(INT x,INT y,INT z,INT a){
		return a+SUNcAlgebra::NumberOfGenerators*Index3D(x,y,z);
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	//  rho(x,y,z) **Not to be confused with overall charge sampling                //
	//  rho(k)=sqrt(CF(k))*Noise(k)                                                 //
	/////////////////////////////////////////////////////////////////////////////////
    
    DOUBLE GammaK(DOUBLE x,DOUBLE kTSqr){
        
        DOUBLE Qs=(Parameters::Qs0A)*pow(x,-0.1)*(1-x)/(sqrt(Lattice::a[0]*Lattice::a[1]));
        
			return 4.0*M_PI*Nc/(Nc*Nc-1)*kTSqr*kTSqr/(Qs*Qs)*exp(-kTSqr/(Qs*Qs));                   
    }
	
	void ChargeFluctuation(ChargeDensity *rho,char D){
		
		//Sample noise \zeta(x)=1/a^3/2*RNG //                                                                  D=[a^(-3/2)]
		FourierVariables= new FFT3D(Lattice::N[0],Lattice::N[1],Lattice::N[2],SUNcAlgebra::NumberOfGenerators);
		
		for(INT z=0;z<=Lattice::N[2]-1;z++){
			for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
					for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						
						DOUBLE StochasticNoise=RandomNumberGenerator::Gauss()/sqrt(Lattice::aCube);
						FourierVariables->SetX(x,y,z,a,StochasticNoise);
					}
				}
			}
		}
		
		
		// PERFORM FFT X->P //
		FourierVariables->ExecuteXtoP();
		
		DOUBLE NormalizationFactor=(1.0)/(Lattice::N[0]*Lattice::N[1]*Lattice::N[2]);
		
		// CF(k)=(N_c/(N_c^2-1))*4\pi*kt^4/Qs^2*exp(-kt^2/Qs^2) //                                             D=[a^-2]
		
		// \rho(x)=\int d^3k sqrt(CF(k))*\zeta(k) exp(ikx) //												   D=[a^-3*a^-1*a^(-3/2)] =[a^-5/2]
		
		for(INT pZIndex=0;pZIndex<=Lattice::N[2]-1;pZIndex++){
			for(INT pYIndex=0;pYIndex<=Lattice::N[1]-1;pYIndex++){
				for(INT pXIndex=0;pXIndex<=Lattice::N[0]-1;pXIndex++){
					
					DOUBLE pLongAbs=(1/M_SQRT2)*sqrt((DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pZIndex/DOUBLE(Lattice::N[2])))/SQR(Parameters::aZInFm));
					
					DOUBLE pTSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pXIndex/DOUBLE(Lattice::N[0])))/SQR(Lattice::a[0])+(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pYIndex/DOUBLE(Lattice::N[1])))/SQR(Lattice::a[1]);
					
					DOUBLE x=pLongAbs/(Parameters::CoMEnergy/M_SQRT2);                                                   //D=[fm^-1/fm^-1]

					for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                        
                        COMPLEX Xi;
						if(x<=1.0){
							
							Xi=NormalizationFactor*FourierVariables->GetP(pXIndex,pYIndex,pZIndex,a)*sqrt(GammaK(x,pTSqr));
							
						}
						else{
							Xi=COMPLEX(0.0,0.0);
						}
                        
                        FourierVariables->SetP(pXIndex,pYIndex,pZIndex,a,Xi);

					}
				}
			}
		}
		
		FourierVariables->ExecutePtoX();
		
		for(INT z=0;z<=Lattice::N[2]-1;z++){
			for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
					for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						
						rho->Get(x,y,z,a)[0]=real(FourierVariables->GetX(x,y,z,a));
					}
				}
			}
		}
		
		delete FourierVariables;
		
	}
}

#endif


