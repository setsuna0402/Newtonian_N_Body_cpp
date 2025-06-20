/***********************************************************************
/
/  Compute acceleration of each particle by soften gravity (CPU version)
/  Note that Gravitational constant is always assumed to be 1 !
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Size (const uint) : Number of particles
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Update position by a given timestep
void UpdatePosition_CPU(real dt, real (*Pos)[3], real (*Vel)[3], const uint Size){
    static uint call_counter = 0;
    if (call_counter == 0) cout << "Which function is running?   " << __func__ << endl;

    // In this case, for-loop parallelisation can be not beneficial. 
    // (Too small workload per thread)
    // However, if you want to use OpenMP, you can uncomment the following lines.
#ifdef OPEN_MP
    #pragma omp parallel
#endif
    {  
    #ifdef OPEN_MP
        #pragma omp for
    #endif
        // Update the position of each particle
        for (uint i = 0; i < Size; i++){
            Pos[i][0] += Vel[i][0] * dt;
            Pos[i][1] += Vel[i][1] * dt;
            Pos[i][2] += Vel[i][2] * dt;
        }
    }
    if (call_counter == 0) cout << __func__ << "...done!" << endl;
    call_counter++;
} // FUNCTION : UpdatePosition_CPU