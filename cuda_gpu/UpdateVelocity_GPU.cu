/***********************************************************************
/
/  Update the velocity of all the particles (GPU version)
/
/  written by: KH 
/  date:       2025/6/22
/  Input:
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Acc  (Dynamical array) : Acceleration of each particle (Size: N x 3)
/  Size (const uint) : Number of particles
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

// GPU function : Update the Velocity of all the particles by a given timestep
__global__ void UpdateVelocity_GPU(real dt, real (*Vel)[3], real (*Acc)[3],
                                   const uint Size){
    // get the target particle index of each thread
    const int p = blockDim.x*blockIdx.x + threadIdx.x;
    if (p < Size) {
        Vel[p][0] += Acc[p][0] * dt;
        Vel[p][1] += Acc[p][1] * dt;
        Vel[p][2] += Acc[p][2] * dt;
    }

} // FUNCTION : UpdateVelocity_GPU