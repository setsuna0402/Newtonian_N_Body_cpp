/***********************************************************************
/
/  Update the position of all the particles (GPU version)
/
/  written by: KH 
/  date:       2025/6/22
/  Input:
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Size (const uint) : Number of particles
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

// GPU function : Update the position of all the particles by a given timestep
__global__ void UpdatePosition_GPU(real dt, real (*Pos)[3], real (*Vel)[3],
                                   const uint Size){
    // get the target particle index of each thread
    const int p = blockDim.x*blockIdx.x + threadIdx.x;
    if (p < Size) {
        Pos[p][0] += Vel[p][0] * dt;
        Pos[p][1] += Vel[p][1] * dt;
        Pos[p][2] += Vel[p][2] * dt;
    }
} // FUNCTION : UpdatePosition_GPU