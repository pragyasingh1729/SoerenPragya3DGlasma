///////////////////////////////////////////////////////////////////
//   IMPLEMENTATION OF ONE DIMENSIONAL FAST FOURIER TRANSFORM     //
///////////////////////////////////////////////////////////////////

#ifndef __FFT_1D_CPP__
#define __FFT_1D_CPP__

//////////////////////////////////////////////
//   TO DO:SET THE PARAMETERS ACCORDINGLY   //
//////////////////////////////////////////////


class FFT1D{
    
    // FFTW PARAMETERS //
private:
    
    static const INT MyDimension=1;
    
    INT MyNx;
    
    INT *MyNumberOfSites;
    
    INT MyNumberOfTransforms;
    
    INT MyOverallSize;
    
    INT MyInputMemDistance;
    INT MyOutputMemDistance;
    
    INT MyInputDataSpacing;
    INT MyOutputDataSpacing;
    
    INT *MyInEmbed;
    INT *MyOutEmbed;
    
    
    // DATA STORAGE //
private:
    fftw_complex *XArray;
    fftw_complex *PArray;
    
    // INDEXING IN COORDINATE AND MOMENTUM SPACE //
    
private:
    INT XIndex(INT x,INT a){
	   return a+MyNumberOfTransforms*(x);
    }
    
    INT PIndex(INT px,INT a){
	   return a+MyNumberOfTransforms*(px);
    }
    
public:
    
    void SetX(INT x,INT a,DOUBLE Value){
	   XArray[XIndex(x,a)][0]=double(Value); XArray[XIndex(x,a)][1]=double(0.0);
    }
    
    void SetP(INT px,INT a,COMPLEX Value){
	   PArray[PIndex(px,a)][0]=double(real(Value)); PArray[PIndex(px,a)][1]=double(imag(Value));
    }
    
    COMPLEX GetX(INT x,INT a){
	   return COMPLEX(XArray[XIndex(x,a)][0],XArray[XIndex(x,a)][1]);
    }
    
    COMPLEX GetP(INT px,INT a){
	   return COMPLEX(PArray[PIndex(px,a)][0],PArray[PIndex(px,a)][1]);
    }
    
    // FFTW PLANS //
private:
    fftw_plan XtoP;
    fftw_plan PtoX;
    
    
    
   
public:
    
    FFT1D(INT Nx,INT NumberOfTransforms){
	   
	   ////////////////////////
	   //SET FFTW PARAMETERS //
	   ////////////////////////
	   
	   //NUMBER OF SITES IN EACH DIMENSION
	   MyNx=Nx;
	   
	   // SET IN REVERSED ORDERING  TO SWITCH FROM ROW-MAJOR TO COLUMN-MAJOR FORMAT //
	   MyNumberOfSites=new INT[1];
	   
	   MyNumberOfSites[0]=MyNx;
	   
	   
	   //NUMBER OF INDEPENDENT TRANSFORMS
	   MyNumberOfTransforms=NumberOfTransforms;
	   
	   //OVERALL DATA SIZE
	   MyOverallSize=MyNx*NumberOfTransforms;
	   
	   //SPACING BETWEEN DATASETS
	   
	   MyInputMemDistance=1; MyOutputMemDistance=1;
	   
	   //SPACING BETWEEN POINTS OF THE SAME DATASET
	   MyInputDataSpacing=MyNumberOfTransforms;  MyOutputDataSpacing=MyNumberOfTransforms;
	   
	   // EMBEDDING OF INDIVIDUAL DATASETS //
	   MyInEmbed=MyNumberOfSites;
	   MyOutEmbed=MyNumberOfSites;
	   
	   ////////////////////////
	   //  ALLOCATE ARRAYS   //
	   ////////////////////////
	   
	   XArray=fftw_alloc_complex(MyOverallSize);
	   PArray=fftw_alloc_complex(MyOverallSize);
	   
	   ////////////////////////
	   //  COMPUTE PLANS     //
	   ////////////////////////
	   
	   
	   XtoP=fftw_plan_many_dft(MyDimension,MyNumberOfSites,MyNumberOfTransforms,
						  XArray,MyInEmbed,
						  MyInputDataSpacing, MyInputMemDistance,
						  PArray,MyOutEmbed,
						  MyOutputDataSpacing,MyOutputMemDistance,
						  FFTW_FORWARD,MY_FFTW_PLANNER_FLAG);
	   
	   PtoX=fftw_plan_many_dft(MyDimension,MyNumberOfSites,MyNumberOfTransforms,
						  PArray,MyOutEmbed,
						  MyOutputDataSpacing, MyOutputMemDistance,
						  XArray,MyInEmbed,
						  MyInputDataSpacing, MyInputMemDistance,
						  FFTW_BACKWARD,MY_FFTW_PLANNER_FLAG);
	   
	   
    }
    
    ///////////////////////
    // CLASS DESTRUCTOR  //
    ///////////////////////
    
    ~FFT1D(){
	   
	   //DESTROY PLANS//
	   fftw_destroy_plan(XtoP);
	   fftw_destroy_plan(PtoX);
	   
	   //DELETE FFTW PARAMETERS //
	   delete[] MyNumberOfSites;
	   
	   //FREE MEMORY //
	   fftw_free(XArray);
	   fftw_free(PArray);
	   
    }
    
    
    //////////////////////////
    //  EXECUTION COMMANDS  //
    //////////////////////////
public:
    
    void ExecuteXtoP(){
	   fftw_execute(XtoP);
    }
    
    void ExecutePtoX(){
	   fftw_execute(PtoX);
    }
    
public:
    
    void OutputX(std::string fname){
	   
	   std::ofstream OutStream;
	   
	   OutStream.open(fname.c_str());
	  
		  for(INT x=0;x<MyNx;x++){
				
				
				OutStream << x << " ";
				
				for(INT a=0;a<MyNumberOfTransforms;a++){
				    
				    OutStream << real(GetX(x,a)) << " " << imag(GetX(x,a)) << " ";
				}
				
				OutStream << std::endl;

	     }
	
	   OutStream.close();
	   
    }
    
    void OutputP(std::string fname){
	   
	   std::ofstream OutStream;
	   
	   OutStream.open(fname.c_str());
	   
		  for(INT px=0;px<MyNx;px++){
		
			 OutStream << px <<  " ";
			 
			 for(INT a=0;a<MyNumberOfTransforms;a++){
				
				OutStream << real(GetP(px,a)) << " " << imag(GetP(px,a)) << " ";
			 }
			 
			 OutStream << std::endl;
			 
	   }
	   
	   OutStream.close();
	   
    }
    
};

#endif
