#ifndef __SUNCVARIABLES__CPP__
#define __SUNCVARIABLES__CPP__

namespace SUNcGaugeLinks{
	
    GaugeLinks *U,*UOld;
    
	void Init(){

	   UOld=new GaugeLinks(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   U=new GaugeLinks(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
    }
    
}

namespace SUNcWilsonLines{
    
    WilsonLines *VRight,*VLeft;
    
    void Init(){
		
	   VRight=new WilsonLines(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   VLeft=new WilsonLines(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
    }
}

namespace SUNcElectricFields{
    
    ElectricFields *E, *B;
    
    void Init(){
	   
	   E=new ElectricFields(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   B=new ElectricFields(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
    }
    
}
 
namespace SUNcChargeDensity {

    ChargeDensity *Jp, *Jm, *JpStat, *JmStat, *rhoL,* rhoR, *RhoRT1, *RhoLT1,*JmForward,*rhoOldL,*rhoOldR,*CFL, *CFR;

    void Init(){

	   rhoL=new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   rhoR=new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   Jp=new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   Jm=new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   JmStat=new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   JpStat=new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   RhoRT1 =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   RhoLT1 =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   JmForward =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   rhoOldL =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   rhoOldR =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   CFL =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
	   CFR =new ChargeDensity(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);

	   
    }

}


namespace Lattice{
    
    ///////////////////////////
    //INITIALIZATION ROUTINE //
    ///////////////////////////
    
    void Init(){
	   
	   // COMMANDLINE OUTPUT OF BOUNDARY CONDITIONS //
	   #if (BOUNDARY_CONDITIONS_FLAG==OPEN_BOUNDARY_CONDITIONS_FLAG)
	   std::cerr << "#SETTING UP LATTICE WITH OPEN BOUNDARY CONDITIONS" << std::endl;
	   #endif
	   
        #if (BOUNDARY_CONDITIONS_FLAG==PERIODIC_BOUNDARY_CONDITIONS_FLAG)
	   std::cerr << "#SETTING UP LATTICE WITH PERIODIC BOUNDARY CONDITIONS" << std::endl;
        #endif
        
        Volume=Lattice::N[0]*Lattice::N[1]*Lattice::N[2];
        aCube =Lattice::a[0]*Lattice::a[1]*Lattice::a[2];
        
        SUNcGaugeLinks::Init();
        SUNcWilsonLines::Init();
        SUNcElectricFields::Init();
        SUNcChargeDensity::Init();
        
    }
    
}




#endif
