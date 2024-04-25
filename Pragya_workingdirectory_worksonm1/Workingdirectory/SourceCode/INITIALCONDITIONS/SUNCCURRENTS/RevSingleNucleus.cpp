namespace InitialConditions{
    
    FFT1D *FourierVariable;
    
    namespace SingleSheet{
	   
	   void SingleSheet(ChargeDensity *rho, INT zLow,INT zHigh,ChargeDensity *RhoEvolved,WilsonLines *V,char D,ChargeDensity *JStat,ChargeDensity *J){

	   rho->SynchronizeGhostCells();
	   V->SynchronizeGhostCells();

		  
		 // GAUGE LINKS AT TIME t //
		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
				    
				    //COMPUTE MATRIX PRODUCTS
				    SUNcGroup::Operations::DU(V->Get(x,y,z),V->Get(x+1,y,z),SUNcGaugeLinks::UOld->Get(x,y,z,0));
				    SUNcGroup::Operations::DU(V->Get(x,y,z),V->Get(x,y+1,z),SUNcGaugeLinks::UOld->Get(x,y,z,1));
				    
				    COPY_SUNcMatrix(SUNcGaugeLinks::UOld->Get(x,y,z,2),SUNcGroup::UnitMatrix);
				}
			 }
		  }
		  
		  SUNcGaugeLinks::UOld->SynchronizeGhostCells(D);
		  
		  //////////////////////////////////////////////////////////////////////////////////
		  //   Rho  EVOLVED USING FOURIER HELPS IN SETTING THE INITIAL ELECTRIC FIELDS   //
		  /////////////////////////////////////////////////////////////////////////////////
		  
		  Dynamics::SUNc::EvolveRhoViaFourier(0,Lattice::N[0]-1,0,Lattice::N[1]-1,0,Lattice::N[2]-1,RhoEvolved,rho,Dynamics::dTStep,D);
		  
		  RhoEvolved->SynchronizeGhostCells();
		  
		  CollidingObjects::SetWilsonLinesForFirstStep(V,RhoEvolved,D);
		  
		//   TmdBasedCollision::SetWilsonLinesForFirstStep(V,RhoEvolved,D);
		   
		  ///////////////////////////////
		  // GAUGE LINKS AT TIME t+dt //
		  //////////////////////////////
		  
		 V->SynchronizeGhostCells();

		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
				    
				    //COMPUTE MATRIX PRODUCTS
				    SUNcGroup::Operations::DU(V->Get(x,y,z),V->Get(x+1,y,z),SUNcGaugeLinks::U->Get(x,y,z,0));
				    SUNcGroup::Operations::DU(V->Get(x,y,z),V->Get(x,y+1,z),SUNcGaugeLinks::U->Get(x,y,z,1));
				    
				    COPY_SUNcMatrix(SUNcGaugeLinks::U->Get(x,y,z,2),SUNcGroup::UnitMatrix);
				}
			 }
		  }
		  
		  SUNcGaugeLinks::U->SynchronizeGhostCells(D);
		  
		  /////////////////////////////////////////////////////
		  //   INITIALIZE ELECTRIC FIELDS  AT TIME (t+dt/2)  //
		  ////////////////////////////////////////////////////
		  
		  DOUBLE coff[Lattice::Dimension-1];
		  
		  for (INT mu=0;mu<Lattice::Dimension-1;mu++){
			 coff[mu]= -(Dynamics::dTStep)*SQR(SUNcGaugeLinks::U->a[mu])/(SUNcGaugeLinks::U->aCube);
		  }
		  
		  SU_Nc_FUNDAMENTAL_FORMAT XProd[SUNcGroup::MatrixSize];
		  SU_Nc_FUNDAMENTAL_FORMAT YProd[SUNcGroup::MatrixSize];
		  
		  for(INT z=zLow;z<=zHigh;z++){
			 for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
				    
				    SUNcGroup::Operations::UD(SUNcGaugeLinks::U->Get(x,y,z,0),SUNcGaugeLinks::UOld->Get(x,y,z,0),XProd);
				    SUNcGroup::Operations::UD(SUNcGaugeLinks::U->Get(x,y,z,1),SUNcGaugeLinks::UOld->Get(x,y,z,1),YProd);
				    
				    SUNcAlgebra::Operations::MatrixILog(coff[0],XProd,SUNcElectricFields::E->Get(x,y,z,0,0));
				    SUNcAlgebra::Operations::MatrixILog(coff[1],YProd,SUNcElectricFields::E->Get(x,y,z,1,0));
				    
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   SUNcElectricFields::E->Get(x,y,z,2,a)[0]=DOUBLE(0.0);
				    }
				}
			 }
		  }
		  
		  SUNcElectricFields::E->SynchronizeGhostCells(D);
		  
		  /////////////////////////////////////////////////////////////
		  //   CALCULATING JSTAT USING GAUSS LAW AT (t+dt/2)        //
		  ///////////////////////////////////////////////////////////
		  
		  SU_Nc_ALGEBRA_FORMAT ExDown[SUNcAlgebra::NumberOfGenerators];
		  SU_Nc_ALGEBRA_FORMAT EyDown[SUNcAlgebra::NumberOfGenerators];
		  
		  for(INT z=0;z<=Lattice::N[2]-1;z++){
			 for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
				    
				    if(z>=zLow && z<=zHigh){

					   SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcGaugeLinks::U->Get(x-1,y,z,0),SUNcElectricFields::E->Get(x-1,y,z,0,0),ExDown);
					   SUNcAlgebra::Operations::InverseAdjointMultiplication(SUNcGaugeLinks::U->Get(x,y-1,z,1),SUNcElectricFields::E->Get(x,y-1,z,1,0),EyDown);
					   
					   for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						  JStat->Get(x,y,z,a)[0]=M_SQRT2*((SUNcElectricFields::E->GetValue(x,y,z,0,a)-ExDown[a])+(SUNcElectricFields::E->GetValue(x,y,z,1,a)-EyDown[a]));
					   }
				    }
    
				    else{
					   for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
						  JStat->Get(x,y,z,a)[0]=0.0;
					   }
				    }
				    
				}
			 }
		  }
		  
		  JStat->SynchronizeGhostCells();
		  
		  //////////////////////////////////////////////////////////////////
		  //  CALCULATING Js USING FOURIER TRANSFORMATION  AT J(t+dt)	//
		  /////////////////////////////////////////////////////////////////
		  
		  FourierVariable=new FFT1D(Lattice::N[2],SUNcAlgebra::NumberOfGenerators*Lattice::N[0]*Lattice::N[1]);

		  for(INT z=0;z<=Lattice::N[2]-1;z++){
			 for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   
					   FourierVariable->SetX(z,SUNcAlgebraIndex2D(x,y,a),JStat->GetValue(x,y,z,a));
				    }
				}
			 }
		  }
		  
		  FourierVariable->ExecuteXtoP();
		  
		  DOUBLE NormalizationFactor=(1.0)/(Lattice::N[2]);
		  
		  for(INT pZIndex=0;pZIndex<=Lattice::N[2]-1;pZIndex++){
			 
			 DOUBLE OmegaSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2])))/SQR(Lattice::a[2]);
			 COMPLEX pZD;
			 COMPLEX pNpos;
			 COMPLEX pNneg;
			 
		  if(D=='R'){
			 pZD=(DOUBLE(1.0)-std::exp(-2.0*ComplexI*M_PI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2])))/Lattice::a[2];
			 pNpos=(std::exp(+ComplexI*sqrt(OmegaSqr)*Dynamics::dTStep)-DOUBLE(1.0));
			 pNneg=(std::exp(-ComplexI*sqrt(OmegaSqr)*Dynamics::dTStep)-DOUBLE(1.0));
		  }
		  
		  if(D=='L'){
			 pZD=-(DOUBLE(1.0)-std::exp(-2.0*ComplexI*M_PI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2])))/Lattice::a[2];
			 pNpos=(std::exp(-ComplexI*sqrt(OmegaSqr)*Dynamics::dTStep)-DOUBLE(1.0));
			 pNneg=(std::exp(+ComplexI*sqrt(OmegaSqr)*Dynamics::dTStep)-DOUBLE(1.0));
		  }

		  for(INT y=0;y<=Lattice::N[1]-1;y++){
			 for(INT x=0;x<=Lattice::N[0]-1;x++){
				for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
				    
				    COMPLEX JStatTilde=NormalizationFactor*FourierVariable->GetP(pZIndex,SUNcAlgebraIndex2D(x,y,a));
				    
				    COMPLEX JTilde;
				    
				    if(pZIndex==0){
					   JTilde=JStatTilde;
				    }
				    else if(pZIndex<Lattice::N[2]/2){
					   JTilde=(JStatTilde*pNpos)/(Dynamics::dTStep*pZD);
				    }
				    else if(pZIndex==Lattice::N[2]/2){
					   JTilde=COMPLEX(0.0,0.0);
				    }
				    else {
					   JTilde=(JStatTilde*pNneg)/(Dynamics::dTStep*pZD);
				    }
				    
				    FourierVariable->SetP(pZIndex,SUNcAlgebraIndex2D(x,y,a),JTilde);
				    }
				}
			 }
		  }
		  
		  FourierVariable->ExecutePtoX();
		  
		  for(INT z=0;z<=Lattice::N[2]-1;z++){
			 for(INT y=0;y<=Lattice::N[1]-1;y++){
				for(INT x=0;x<=Lattice::N[0]-1;x++){
				    for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
					   
					   J->Get(x,y,z,a)[0]=real(FourierVariable->GetX(z,SUNcAlgebraIndex2D(x,y,a)));
				    }
				}
			 }
		  }
		  
		  delete FourierVariable;
		  
		  J->SynchronizeGhostCells();
	   }
    }
}
		  
