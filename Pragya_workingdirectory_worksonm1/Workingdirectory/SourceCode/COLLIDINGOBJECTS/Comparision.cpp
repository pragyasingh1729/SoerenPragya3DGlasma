namespace ComparisionCheck{
    
    void JmConstructedProfile(){

	   std::ofstream JmLSitOutStream;
	   JmLSitOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"JmCheck",".txt").c_str());

	   // SET PRECISION //
	   JmLSitOutStream.precision(OUTPUT_PRECISION);

	   Dynamics::SUNc::EvolveRhoViaFourier(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,SUNcChargeDensity::JmForward,SUNcChargeDensity::Jm ,0.5*Dynamics::dTStep,'L');

	   // CREATE OUTPUT //
	   for(INT z=0;z<=Lattice::N[2]-1;z++){

		  DOUBLE AvgJmBack[SUNcAlgebra::NumberOfGenerators];
		  SU_Nc_ALGEBRA_FORMAT AvgVJmVD[SUNcAlgebra::NumberOfGenerators];
		  SU_Nc_ALGEBRA_FORMAT VJmVD[SUNcAlgebra::NumberOfGenerators];

		  for (INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
			 AvgVJmVD[a]=0.0;
			 AvgJmBack[a]=0.0;
		  }

		  for(INT y=0;y<=Lattice::N[1]-1;y++){
			 for(INT x=0;x<=Lattice::N[0]-1;x++){

				SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcWilsonLines::VRight->Get(x,y,z),SUNcChargeDensity::RhoRT1->Get(x,y,z,0),VJmVD);

				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    AvgVJmVD[a]+=VJmVD[a];
				    AvgJmBack[a]+=SUNcChargeDensity::JmForward->GetValue(x,y,z,0);
				}
			 }
		  }

		  JmLSitOutStream << z << " " <<  (AvgVJmVD[0])/(Lattice::N[0]*Lattice::N[1]) <<" " << (AvgJmBack[0])/(Lattice::N[0]*Lattice::N[1])<< std::endl;
	   }

	   //CLOSE OUTPUT STREAM
	   JmLSitOutStream.close();
    }


    void JpConstructedProfile(){
	   
	   std::ofstream JpRSitOutStream;
	   JpRSitOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"JpCheck",".txt").c_str());
	   
	   // SET PRECISION //
	   JpRSitOutStream.precision(OUTPUT_PRECISION);
	   
	   // CREATE OUTPUT //
	   for(INT z=0;z<=Lattice::N[2]-1;z++){
		  
		  DOUBLE AvgJpStat=0.0;
		  DOUBLE AvgJp=0.0;

		  SU_Nc_ALGEBRA_FORMAT AvgVJpVD[SUNcAlgebra::NumberOfGenerators];
		  SU_Nc_ALGEBRA_FORMAT VJpVD[SUNcAlgebra::NumberOfGenerators];
		  
		  for (INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
			 AvgVJpVD[a]=0.0;
		  }
		  
		  for(INT y=0;y<=Lattice::N[1]-1;y++){
			 for(INT x=0;x<=Lattice::N[0]-1;x++){
				
				SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcWilsonLines::VLeft->Get(x,y,z),SUNcChargeDensity::RhoLT1->Get(x,y,z,0),VJpVD);
				
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    AvgVJpVD[a]+=VJpVD[a];
				    AvgJpStat+=SUNcChargeDensity::JpStat->Get(x,y,z,0)[0];
				    AvgJp+=SUNcChargeDensity::Jp->Get(x,y,z,0)[0];

				}
			 }
		  }
		  
		  JpRSitOutStream << z << " " << (AvgVJpVD[0])/(Lattice::N[0]*Lattice::N[1])<< " " << (AvgJpStat)/(Lattice::N[0]*Lattice::N[1])<< " " << (AvgJp)/(Lattice::N[0]*Lattice::N[1])<<std::endl;
	   }
	   
	   //CLOSE OUTPUT STREAM
	   JpRSitOutStream.close();
    }


    void JmFullConstructedProfile(){

	   std::ofstream FullJmLSitOutStream; 
	   FullJmLSitOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"JmCheck",".txt").c_str());

	   // SET PRECISION //
	   FullJmLSitOutStream.precision(OUTPUT_PRECISION);

	   // CREATE OUTPUT //
	   Dynamics::SUNc::EvolveRhoViaFourier(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,SUNcChargeDensity::JmForward,SUNcChargeDensity::JmStat,4000*Dynamics::dTStep,'R');

	   SU_Nc_ALGEBRA_FORMAT VJmVD[SUNcAlgebra::NumberOfGenerators];

	   for(INT z=-1;z<=Lattice::N[2];z++){

		  DOUBLE TrJfJf=0.0;   DOUBLE JfJf=0.0;
		  DOUBLE TrJmJm=0.0;   DOUBLE JmJm=0.0;
		  DOUBLE VRhoVDSqr=0;
		   
		   for(INT y=0;y<=Lattice::N[1]-1;y++){
			  for(INT x=0;x<=Lattice::N[0]-1;x++){
				  for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					
					SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcWilsonLines::VRight->Get(x,y,z),SUNcChargeDensity::RhoRT1->Get(x,y,z,0),VJmVD);

				    VRhoVDSqr+=SQR(VJmVD[a])/(Lattice::N[0]*Lattice::N[1]);
					JfJf+=SQR(SUNcChargeDensity::JmForward->Get(x,y,z,a)[0]);
					JmJm+=SQR(SUNcChargeDensity::Jm->Get(x,y,z,a)[0]);

				}
			 }
		  }
		
			TrJmJm=DOUBLE(0.5)*(JmJm)/(Lattice::N[0]*Lattice::N[1]);
			TrJfJf=DOUBLE(0.5)*(JfJf)/(Lattice::N[0]*Lattice::N[1]);
			
		   FullJmLSitOutStream << z << " " << TrJmJm<< " "<< DOUBLE(0.5)*VRhoVDSqr << "  "<< TrJfJf <<std::endl;

	   }

	   //CLOSE OUTPUT STREAM
	   FullJmLSitOutStream.close();
    }
    

   void JmConstructedComp(){

	   std::ofstream JmLSitOutStream;
	   JmLSitOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"JmColorComp",".txt").c_str());

	   // SET PRECISION //
	   JmLSitOutStream.precision(OUTPUT_PRECISION);

	   // CREATE OUTPUT //
	   Dynamics::SUNc::EvolveRhoViaFourier(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,SUNcChargeDensity::JmForward,SUNcChargeDensity::JmStat,4000*Dynamics::dTStep,'R');

	   SU_Nc_ALGEBRA_FORMAT VJmVD[SUNcAlgebra::NumberOfGenerators];

	   for(INT z=-1;z<=Lattice::N[2];z++){
	       DOUBLE JfJf=0.0;  DOUBLE JmJm=0.0;  DOUBLE VRhoVDSqr=0;
			   for(INT y=0;y<=Lattice::N[1]-1;y++){
				  for(INT x=0;x<=Lattice::N[0]-1;x++){
					
					SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcWilsonLines::VRight->Get(x,y,z),SUNcChargeDensity::RhoRT1->Get(x,y,z,0),VJmVD);

					VRhoVDSqr+=SQR(VJmVD[0])/(Lattice::N[0]*Lattice::N[1]);
					JfJf+=SQR(SUNcChargeDensity::JmForward->GetValue(x,y,z,0))/(Lattice::N[0]*Lattice::N[1]);
					JmJm+=SQR(SUNcChargeDensity::Jm->GetValue(x,y,z,0))/(Lattice::N[0]*Lattice::N[1]);

				}
			 }
			
			JmLSitOutStream << z << " " << JmJm<< " "<< VRhoVDSqr << "  "<< JfJf <<std::endl;

	   }

	   //CLOSE OUTPUT STREAM
	   JmLSitOutStream.close();
	}


    // CHECK [Di,Fiz]=Jz //
    void YMECheckforJz(){
	   std::ofstream YMECheckOutStream;
	   YMECheckOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"YME_JzCheck.txt").c_str());
	   
	   // SET PRECISION //
	   YMECheckOutStream.precision(OUTPUT_PRECISION);
	   
	   DOUBLE DiFiz[SUNcAlgebra::NumberOfGenerators];
	   SU_Nc_ALGEBRA_FORMAT BxDown[SUNcAlgebra::NumberOfGenerators];
	   SU_Nc_ALGEBRA_FORMAT ByDown[SUNcAlgebra::NumberOfGenerators];
	 
	   for(INT z=0;z<=Lattice::N[2]-1;z++){
		  for(INT y=0;y<=Lattice::N[1]-1;y++){
			 for(INT x=0;x<=Lattice::N[0]-1;x++){
				
				SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcGaugeLinks::U->Get(x,y-1,z,1),SUNcElectricFields::B->Get(x,y-1,z,0,0),BxDown);
				SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcGaugeLinks::U->Get(x-1,y,z,0),SUNcElectricFields::B->Get(x-1,y,z,1,0),ByDown);
				
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    DiFiz[a]= -(SUNcElectricFields::B->GetValue(x,y,z,1,a)-ByDown[a])+(SUNcElectricFields::B->GetValue(x,y,z,0,a)-BxDown[a]);
				}
				
				YMECheckOutStream << x << " " << y << " " << z << " " << DiFiz[0] <<"  " << (SUNcChargeDensity::Jp->GetValue(x,y,z,0)-SUNcChargeDensity::Jm->GetValue(x,y,z,0))/M_SQRT2 << std::endl;

			 }
		  }
	   }
	   
    }
}



