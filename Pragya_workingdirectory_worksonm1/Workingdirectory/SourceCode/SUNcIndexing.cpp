#ifndef __SUNC_INDEXING__CPP__
#define __SUNC_INDEXING__CPP__
    
//GENERAL 2D INDEX
#define Index2D(x,y)  (MOD((x),Lattice::N[0])+Lattice::N[0]*(MOD((y),Lattice::N[1])))

#define Index3D(x,y,z) (MOD(x,Lattice::N[0])+Lattice::N[0]*(MOD(y,Lattice::N[1])+Lattice::N[1]*(z+1)))

#define VectorFieldIndex(x,y,z,mu,a) ((a)+SUNcAlgebra::NumberOfGenerators*((mu)+Lattice::Dimension*Index3D((x),(y),(z))))

#define VectorIndex(x,y,z,a) ((a)+SUNcAlgebra::NumberOfGenerators*(Index3D((x),(y),(z))))

#define GaugeLinkIndex(x,y,z,mu)  (SUNcGroup::MatrixSize*((mu)+Lattice::Dimension*Index3D((x),(y),(z))))


//SU(Nc) MATRIX INDEXING
    INT SUNcMatrixIndex2D(INT x,INT y){
        return SUNcGroup::MatrixSize*(Index2D(x,y));
    }

    
    //SU(Nc) ALGEBRA INDEXING
    INT SUNcAlgebraIndex2D(INT x,INT y,INT a){
        return a+SUNcAlgebra::NumberOfGenerators*(Index2D(x,y));
    }

    
    //TRANSVERSE GAUGE LINK INDEXING
    INT GaugeLinkIndex2D(INT x,INT y,INT mu){
        return SUNcGroup::MatrixSize*(mu+2*(Index2D(x,y)));
    }

#endif


