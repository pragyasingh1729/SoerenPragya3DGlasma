# SoerenPragya3DGlasma


This code was majorly written by Soeren for 2D glasma, I added the 3D glasma part to it. In the current stage, we run three different simulations: collision, left moving, and right moving, which leads to 3 different main files. Atm, one has to manually comment the other two main files depending on the simulation type. 

Apart from that, this code works well on MAC M1 and M2. I have added a header file called `sse2neon.h` which is a translator of Intel SSE (Streaming SIMD Extensions) intrinsics to Arm NEON. It is a replacement of #include <xmmintrin.h>  or #include <emmintrin.h> libraries. So I have replaced them with the #include "sse2neon.h"



