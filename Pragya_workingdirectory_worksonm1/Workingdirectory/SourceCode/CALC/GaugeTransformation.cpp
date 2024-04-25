#ifndef __CALC__GAUGETRANSFORMATION__CPP__
#define __CALC__GAUGETRANSFORMATION__CPP__

//COMPUTE GAUGE TRANSFORMATION OF GAUGE LINKS
#define PERFORM_GAUGE_LINK_TRANSFORMATION(x,y,z)\
\
SUNcGroup::AdvancedOperations::UUD(G->Get(x,y,z),UOld->Get(x,y,z,0),G->Get(x+1,y,z),UNew->Get(x,y,z,0)); \
SUNcGroup::AdvancedOperations::UUD(G->Get(x,y,z),UOld->Get(x,y,z,1),G->Get(x,y+1,z),UNew->Get(x,y,z,1)); \
SUNcGroup::AdvancedOperations::UUD(G->Get(x,y,z),UOld->Get(x,y,z,2),G->Get(x,y,z+1),UNew->Get(x,y,z,2));


//COMPUTE GAUGE TRANSFORMATION OF ELECTRIC FIELDS
#define PERFORM_ELECTRIC_FIELD_TRANSFORMATION(x,y,z) \
\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),EOld->Get(x,y,z,0,0),ENew->Get(x,y,z,0,0));\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),EOld->Get(x,y,z,1,0),ENew->Get(x,y,z,1,0));\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),EOld->Get(x,y,z,2,0),ENew->Get(x,y,z,2,0));


//COMPUTE GAUGE TRANSFORMATION OF COLOR CURRENTS
#define PERFORM_COLOR_CURRENT_TRANSFORMATION(x,y,z) \
\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),jOld->Get(x,y,z,0,0),jNew->Get(x,y,z,0,0));\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),jOld->Get(x,y,z,1,0),jNew->Get(x,y,z,1,0));\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),jOld->Get(x,y,z,2,0),jNew->Get(x,y,z,2,0));


//COMPUTE GAUGE TRANSFORMATION OF COLOR CURRENTS
#define PERFORM_COLOR_CHARGE_TRANSFORMATION(x,y,z) \
\
SUNcAlgebra::Operations::AdjointMultiplication(G->Get(x,y,z),j0Old->Get(x,y,z,0),j0New->Get(x,y,z,0));\

#endif

