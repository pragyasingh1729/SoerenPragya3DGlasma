namespace Timing{
    
    DOUBLE c0;
    
    double clock(){
        return omp_get_wtime();
    }
    
    void Reset(){
        c0=clock();
    }

    
    DOUBLE Get(){
        return clock()-c0;
    }
}
