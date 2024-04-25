namespace InitialConditions{
    
    namespace SUNcGaugeFields{
        
        ///////////////////////////////////////////////
        // SET ALL SU(Nc) FIELDS IDENTICALLY TO ZERO //
        ///////////////////////////////////////////////
        
        void SetZero(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //SET ELECTRIC FIELDS TO ZERO
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                                E->Get(x,y,z,mu,a)[0]=0.0;
                            }
                        }
                        
                        //SET GAUGE LINKS TO UNITY
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            COPY_SUNcMatrix(U->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                        }
                        
                    }
                }
            }
        }
        
        void SetZero(GaugeLinks *U,ElectricFields *E){
            
            SetZero(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
            
        }
        
        void SetZero(){
            
            // SET INITIAL CONDITIONS  //
                std::cerr << "#SETTING PURE GAUGE FIELD INITIAL CONDITIONS" << std::endl;
            
            SetZero(SUNcGaugeLinks::U,SUNcElectricFields::E);
            
            // CHECK INITIAL OBSERVABLES //
                std::cerr << "#CHECKING INITIAL STATE OBSERVABLES" << std::endl;
            
            //CHECK GAUSS LAW VIOLATION //
            Observables::GaussLaw::CheckViolation();

            std::cerr << "#CHECKING UNITARITY VIOLATION" << std::endl;
            
            //CHECK UNITARITY VIOLATION //
            Observables::Unitarity::CheckViolation();
            
            //CHECK INITIAL ENERGY DENSITY //
            Observables::Bulk::Update();
            
                std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
                std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
        
        }
        
        
        ////////////////////////////////////////////////////////////////////////
        // SET SU(Nc) ELECTRIC FIELDS TO ZERO AND GAUGE LINKS TO A PURE GAUGE //
        ////////////////////////////////////////////////////////////////////////
        
        void SetPureGauge(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
            
            // INITIALIZE GAUGE LINKS AND ELECTRIC FIELDS //
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        // SET ELECTRIC FIELDS TO ZERO //
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                                SUNcElectricFields::E->Get(x,y,z,mu,a)[0]=0.0;
                            }
                        }
                        
                        // SET GAUGE LINKS TO UNITY //
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            COPY_SUNcMatrix(SUNcGaugeLinks::U->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                        }
                        
                        
                    }
                }
            }
            
            // RANDOM PURE GAUGE TRANSFORM LINKS //
            GaugeTransformation::SetRandom();
            
            GaugeTransformation::Operations::GaugeTransformLinks();
            GaugeTransformation::Operations::GaugeTransformElectricFields();
            
            Copy(SUNcGaugeLinks::U,GaugeFixedVariables::U);
            Copy(SUNcElectricFields::E,GaugeFixedVariables::E);
            
            
            
        }
        
        
        void SetPureGauge(){
            
            //   SET INITIAL CONDITIONS  //
                std::cerr << "#SETTING PURE GAUGE FIELD INITIAL CONDITIONS" << std::endl;
            
            SetPureGauge(0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1);
            
            // CHECK INITIAL OBSERVABLES //
                std::cerr << "#CHECKING INITIAL STATE OBSERVABLES" << std::endl;
            
            //CHECK GAUSS LAW VIOLATION //
            Observables::GaussLaw::CheckViolation();
            
            //CHECK UNITARITY VIOLATION //
            Observables::Unitarity::CheckViolation();
            
            //CHECK INITIAL ENERGY DENSITY //
            Observables::Bulk::Update();
            
            
                std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
                std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
            
            
        }
        
        
        //////////////////////////////////////////////////////////////////////////////////
        // SET SU(Nc) ELECTRIC FIELDS TO ZERO AND GAUGE LINKS TO RANDOM SU(Nc) MATRICES //
        //////////////////////////////////////////////////////////////////////////////////
        
        
        void SetRandomMatrices(DOUBLE Amplitude,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //SET ELECTRIC FIELDS TO ZERO
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                                SUNcElectricFields::E->Get(x,y,z,mu,a)[0]=0.0;
                            }
                        }
                        
                        //SET GAUGE LINKS TO RANDOM MATRICES
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            RandomNumberGenerator::SUNcMatrix(Amplitude,SUNcGaugeLinks::U->Get(x,y,z,mu));
                        }
                        
                    }
                }
            }
        }
        
        
        void SetRandomMatrices(DOUBLE Amplitude){
            
            //   SET INITIAL CONDITIONS  //
            std::cerr << "#SETTING RANDOM MATRIX INITIAL CONDITIONS" << std::endl;
            
            SetRandomMatrices(Amplitude,0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1);
            
            // CHECK INITIAL OBSERVABLES //
            std::cerr << "#CHECKING INITIAL STATE OBSERVABLES" << std::endl;
        
            //CHECK GAUSS LAW VIOLATION //
            Observables::GaussLaw::CheckViolation();
            
            //CHECK UNITARITY VIOLATION //
            Observables::Unitarity::CheckViolation();
            
            //CHECK INITIAL ENERGY DENSITY //
            Observables::Bulk::Update();
            
                std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
                std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
        
            
        }
        
        void SetRandomMatrices(){
            SetRandomMatrices(4.0*M_PI);
        }
        
        
        ///////////////////////////////////////////////////////////////////
        // SET CONSTANT SU(Nc) ELECTRIC FIELDS AND GAUGE LINKS TO UNITY  //
        ///////////////////////////////////////////////////////////////////
        
        void SetConstantElectricField(DOUBLE Ex,DOUBLE Ey,DOUBLE Ez,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
            
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //SET ELECTRIC FIELDS TO ZERO
                        for(int a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            
                            E->Get(x,y,z,0,a)[0]=Ex;
                            E->Get(x,y,z,1,a)[0]=Ey;
                            E->Get(x,y,z,2,a)[0]=Ez;
                            
                        }
                        
                        //SET GAUGE LINKS TO UNITY
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            COPY_SUNcMatrix(U->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                        }
                        
                    }
                }
            }
            
        }
        
        void SetConstantElectricField(DOUBLE Ex,DOUBLE Ey,DOUBLE Ez,GaugeLinks *U,ElectricFields *E){
            
            SetConstantElectricField(Ex,Ey,Ez,0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
            
        }
        
        void SetConstantElectricField(DOUBLE Ex,DOUBLE Ey,DOUBLE Ez){
            SetConstantElectricField(Ex,Ey,Ez,SUNcGaugeLinks::U,SUNcElectricFields::E);
        }
        
        
        
        
    }
    
    
    
}
