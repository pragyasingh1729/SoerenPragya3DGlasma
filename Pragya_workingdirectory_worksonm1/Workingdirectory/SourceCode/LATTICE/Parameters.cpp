#ifndef __PARAMETERS__CPP__
#define __PARAMETERS__CPP__

namespace Parameters{
	
	INT LeftNucleus=1;		INT RightNucleus=1;

	DOUBLE Qs0A=0.33/sqrt(Lattice::a[0]*Lattice::a[1]);

	//********** with different value of Qs0A for right and left nuclei ****//
	DOUBLE Qs0A_L=0.125/sqrt(Lattice::a[0]*Lattice::a[1]);

	DOUBLE Qs0A_R=0.0625/sqrt(Lattice::a[0]*Lattice::a[1]);

	//********** with different value of Qs0A for right and left nuclei ****//

	DOUBLE Qs0R=1.0;

    // UV REGULARISATIONS //
    DOUBLE Lambda=5*Qs0A;
    
    // IR CUT-OFF //
	DOUBLE MEff=Qs0A;
    
    DOUBLE SigmaLOverAzL=16;
	
	DOUBLE SigmaLOverAzR=16;
	
	//EXPLICITLY FOR PHYSICAL CHARGES//
	//Transverse scales are Qs,a_\perp and R_p. With R_p=2GeV^-1, two dimensionless combination of scales are QsR_p and Qsa_\perp//
	
	DOUBLE SigmaLOverAz=512;
	
	DOUBLE Qs0=0.5; //in Gev calucated from QsRp=1.0 //
			
	INT AtomicNumber=197;
	
	DOUBLE CoMEnergy=500/0.197;												 //in fm^-1
		
	DOUBLE aXInFm=Qs0A*0.197/Qs0;
	
	DOUBLE aYInFm=Qs0A*0.197/Qs0;
	
	DOUBLE aZInFm=aXInFm*Qs0R/(Qs0A*SigmaLOverAz);
	
	DOUBLE aCubeInFm=aXInFm*aYInFm*aZInFm;

	INT Load=0;
	
}

#endif

//** old simulations done for QsA=1/8  R/a_||=16  N_\perp =128 //

