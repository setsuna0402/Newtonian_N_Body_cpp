/***********************************************************************
/
/  Update the position of all the particles (CPU version)
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Acc  (Dynamical array) : Acceleration of each particle (Size: N x 3)
/  Size (const uint) : Number of particles
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Update velocity by a given timestep
void UpdateVelocity_CPU(real dt, real (*Vel)[3], real (*Acc)[3], const uint Size){
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
            Vel[i][0] += Acc[i][0] * dt;
            Vel[i][1] += Acc[i][1] * dt;
            Vel[i][2] += Acc[i][2] * dt;
        }
    }

    if (call_counter == 0) cout << __func__ << "...done!" << endl;
    call_counter++;
} // FUNCTION : UpdateVelocity_CPU