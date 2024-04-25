#ifndef __FOURIER_ACCELERATION_GAUGEFIXING__CPP__
#define __FOURIER_ACCELERATION_GAUGEFIXING__CPP__

namespace FourierAccelerationAlgorithm{
    
    DOUBLE alphaStepSize=2.0;
    
    //DOUBLE alphaStepSize=0.3;
    
    /*MONITORING VARIABLES*/
    DOUBLE MaxStepSize,SqrSumStepSize; 
    DOUBLE MaxDeviation,SqrSumDeviation;  
    DOUBLE MinGaugeFunctional,SumGaugeFunctional;
    
    void UpdateGaugeTransformation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,INT pxLow,INT pxHigh,INT pyLow,INT pyHigh,INT pzLow,INT pzHigh,GaugeLinks *UGaugeTransformed,GaugeTransformations *G){
        
        //CONVERGENCE MONITORING         
        MaxDeviation=0.0;   MinGaugeFunctional=(2.0*Nc)*(1.0/SQR(SUNcGaugeLinks::U->a[0])+1.0/SQR(SUNcGaugeLinks::U->a[1])+1.0/SQR(SUNcGaugeLinks::U->a[2])); 
        
        SqrSumDeviation=0.0;    SumGaugeFunctional=0.0;
        
        #pragma omp parallel 
        {
            
            //SET BUFFERS
            SET_GAUGE_DEVIATION_BUFFERS();
            
            // COMPUTE LOCAL DEVIATION //
            #pragma omp for reduction( + : SqrSumDeviation,SumGaugeFunctional) reduction ( max : MaxDeviation) reduction ( min : MinGaugeFunctional) 
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //GET LOCAL GAUGE LINKS
                        GET_LOCAL_GAUGE_VIOLATION(x,y,z);
                        
                        //MONITOR THE QUALITY OF THE GAUGE FIXING 
                        MONITOR_GAUGE_DEVIATION();
                        
                        for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            FourierSpace::GaugeUpdate->SetX(x,y,z,a,DeviationBuffer[a]);
                        }                    
                        
                    }
                    
                }
            }
            
            
        } // END PARALLEL
        
        // COMPUTE FFT // 
        FourierSpace::GaugeUpdate->ExecuteXtoP();
        
        #pragma omp parallel
        {
            DOUBLE pSqr; COMPLEX NewValue; DOUBLE NormalizationFactor=1.0/(SUNcGaugeLinks::U->N[0]*SUNcGaugeLinks::U->N[1]*SUNcGaugeLinks::U->N[2]);
            
            // PERFORM FOURIER ACCELERATION //
            #pragma omp for
            for(INT pZIndex=pzLow;pZIndex<=pzHigh;pZIndex++){
                for(INT pYIndex=pyLow;pYIndex<=pyHigh;pYIndex++){
                    for(INT pXIndex=pxLow;pXIndex<=pxHigh;pXIndex++){
                        
                        /// DETERMINE MOMENTUM //
                        pSqr=(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pXIndex/DOUBLE(SUNcGaugeLinks::U->N[0])))/SQR(SUNcGaugeLinks::U->a[0])+(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pYIndex/DOUBLE(SUNcGaugeLinks::U->N[1])))/SQR(SUNcGaugeLinks::U->a[1])+(DOUBLE(2.0)-DOUBLE(2.0)*cos(DOUBLE(2.0)*M_PI*pZIndex/DOUBLE(SUNcGaugeLinks::U->N[2])))/SQR(SUNcGaugeLinks::U->a[2]);
                        
                        // SET ACCELERATED VERSION //
                        if(pSqr>0.0){
                            
                            for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                                
                                NewValue=NormalizationFactor*FourierSpace::GaugeUpdate->GetP(pXIndex,pYIndex,pZIndex,a)/pSqr;
                                
                                FourierSpace::GaugeUpdate->SetP(pXIndex,pYIndex,pZIndex,a,NewValue);
                            }
                        }
                    }
                }
            }
            
        } // END PARALLEL
        
        // COMPUTE INVERSE FFT //
        FourierSpace::GaugeUpdate->ExecutePtoX();
        
        //CONVERGENCE MONITORING 
        MaxStepSize=0.0; SqrSumStepSize=0.0; 
        
        // SET BUFFERS //
        #pragma omp parallel
        {
            
            SET_GAUGE_UPDATE_BUFFERS();
            
            DOUBLE UpdateBuffer[SUNcAlgebra::NumberOfGenerators];
            
            // PERFORM COORDINATE SPACE UPDATE //
            #pragma omp for reduction( + : SqrSumStepSize)  reduction ( max : MaxStepSize)

            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        // GET UPDATE //
                        for(INT a=0;a<SUNcAlgebra::NumberOfGenerators;a++){
                            UpdateBuffer[a]=real(FourierSpace::GaugeUpdate->GetX(x,y,z,a));
                        }
                        
                        // COMPUTE GAUGE UPDATE //
                        SUNcAlgebra::Operations::MatrixIExp(-alphaStepSize,UpdateBuffer,wTilde);
                        
                        //MONITOR GAUGE UPDATE
                        MONITOR_GAUGE_UPDATE();
                        
                        //UPDATE GAUGE TRANSFORMATION 
                        SUNcGroup::Operations::DU(wTilde,G->Get(x,y,z),gNew);
                        
                        //WRITE TO GLOBAL ARRAY
                        COPY_SUNcMatrix(G->Get(x,y,z),gNew);
                        
                    }
                }
            }
            
        } // END PARALLEL
        
        
    }

    void UpdateGaugeTransformation(GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
        
        //UPDATE EACH CHECKERBOARD PARITY SEPARATELY
        UpdateGaugeTransformation(0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1,0,SUNcGaugeLinks::U->N[0]-1,0,SUNcGaugeLinks::U->N[1]-1,0,SUNcGaugeLinks::U->N[2]-1,UNew,G);
        GaugeTransformation::Operations::GaugeTransformLinks(UOld,UNew,G);
        
        //SET GLOBAL MONITORING VARIABLES
        GlobalMaxStepSize=MaxStepSize;  GlobalMaxDeviation=MaxDeviation; GlobalMinGaugeFunctional=MinGaugeFunctional;
        
        GlobalAvgStepSize=sqrt(SqrSumStepSize/(SUNcGaugeLinks::U->Volume*SUNcAlgebra::NumberOfGenerators));   GlobalAvgDeviation=sqrt(SqrSumDeviation/(SUNcGaugeLinks::U->Volume*SUNcAlgebra::NumberOfGenerators));     GlobalSumGaugeFunctional=SumGaugeFunctional/(SUNcGaugeLinks::U->Volume);
        
    }
    
    void UpdateGaugeTransformation(){
        
        UpdateGaugeTransformation(SUNcGaugeLinks::U,GaugeFixedVariables::U,GaugeTransformation::G);
        
    }
    
}


#endif
