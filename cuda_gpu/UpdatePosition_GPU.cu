/***********************************************************************
/
/  Update the position of all the particles (GPU version)
/
/  written by: KH 
/  date:       2025/6/22
/  Input:
/  dt  (const real) : Timestep for the update
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
    const int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i < Size) {
        Pos[i][0] += Vel[i][0] * dt;
        Pos[i][1] += Vel[i][1] * dt;
        Pos[i][2] += Vel[i][2] * dt;
    }
} // FUNCTION : UpdatePosition_GPU