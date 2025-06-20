/***********************************************************************
/
/  Compute acceleration of each particle by soften gravity (CPU version)
/  Note that Gravitational constant is always assumed to be 1 !
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
/  Acc  (Dynamical array) : Acceleration of each particle (Size: N x 3)
/  Mass (Dynamical array) : Mass of each particle (Size: N)
/  Size (const uint ) : Number of particles


************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Compute the acceleration for each particles by soften gravity
void Compute_ACC_CPU(real (*Pos)[3], real (*Acc)[3], real *Mass, const uint Size)
{
    static uint call_counter = 0;
    if (call_counter == 0)
        cout << "Which function is running?   " << __func__ << endl;
    const float eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0, 0, 0};  //
    real r = 0.0;
    real r3 = 0.0;

    // calculate the acceleration for the i-th particle from all j particle
    // F_ij : Force on the i-th particle caused by the j particle
#ifdef OPEN_MP
#pragma omp parallel
    {
// dr[3], r and r3 are private variables for each thread
#pragma omp for private(r, r3) firstprivate(dr)
#endif
        for (uint i = 0; i < Size; i++){
            // initialize acc to zero
            for (uint d = 0; d < 3; d++) Acc[i][d] = 0.0;
            // loop over all the particles
            // Note Acc_i_i has to be a zero because there is no self gravity
            // Seperate to two loop to numerically follow Newton's second law
            for (uint j = 0; j < i; j++){
                for (uint d = 0; d < 3; d++) dr[d] = Pos[j][d] - Pos[i][d];
            #ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #endif
                r3 = r * r * r;
                // accumulate the acceleration (acc_ij = Mass_j*r_ij/|r_ij|^3)
                for (uint d = 0; d < 3; d++) Acc[i][d] += Mass[j] * dr[d] / r3;
            }
            for (uint j = i + 1; j < Size; j++){
                // Newton's third law: F_ij = -F_ji
                for (uint d = 0; d < 3; d++) dr[d] = -(Pos[i][d] - Pos[j][d]);
            #ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #endif
                r3 = r * r * r;
                // accumulate the acceleration (acc_ij = Mass_j*r_ij/|r_ij|^3)
                for (uint d = 0; d < 3; d++) Acc[i][d] += Mass[j] * dr[d] / r3;
            }
        }
#ifdef OPEN_MP
    }
#endif

    if (call_counter == 0)
        cout << __func__ << "...done!" << endl;
    call_counter++;
} // FUNCTION : Compute_ACC_CPU