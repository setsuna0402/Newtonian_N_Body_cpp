/***********************************************************************
/
/  Compute Kinetic of each particle  (GPU version)
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Mass (Dynamical array) : Mass of each particle (Size: N)
/  E_Kinetic  (Dynamical array) : Kinetic energy of each particle (Size: N)
/  Size (const uint ) : Number of particles


************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"


// GPU function : Compute the kinetic energy for each particles
__global__ void Compute_KineticEnergy_GPU(real (*Vel)[3], real *Mass,
                                          real *E_Kinetic, const uint Size){
    // get the target particle index of each thread
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    real v2 = 0.0f;
    if (i < Size) {
        // Kinetic energy = 0.5 * m * v^2
        v2 = Vel[i][0] * Vel[i][0] + Vel[i][1] * Vel[i][1] + Vel[i][2] * Vel[i][2];
        E_Kinetic[i] = 0.5 * Mass[i] * v2;
    }

} // Function: Compute_KineticEnergy_GPU