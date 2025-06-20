/***********************************************************************
/
/  Compute potential energy of each particle by soften gravity (CPU version)
/  Note that Gravitational constant is always assumed to be 1 !
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
/  Mass (Dynamical array) : Mass of each particle (Size: N)
/  E_Potential  (Dynamical array) : potential energy of each particle (Size: N)
/  Size (const uint ) : Number of particles


************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Compute the potential energy for each particles
void Compute_PotentialEnergy_CPU(real (*Pos)[3], real *E_Potential,
                                 real *Mass, const uint Size){
    static uint call_counter = 0;
    if (call_counter == 0)
        cout << "Which function is running?   " << __func__ << endl;
    const float eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0, 0, 0};  //
    real r = 0.0;

    // calculate the potential energy for the i-th particle from all j particle
    // dr_ij : Position vector points to j-th particle from the i-th particle
#ifdef OPEN_MP
#pragma omp parallel
    {
// dr[3] and r are private variables for each thread
#pragma omp for private(r) firstprivate(dr)
#endif
        for (uint i = 0; i < Size; i++){
            // initialise PE to zero
            E_Potential[i] = 0.0;
            // loop over all the particles
            // Note E_Potential_i_i has to be a zero because there is no self gravity
            // Seperate to two loop to numerically follow Newton's second law
            for (uint j = 0; j < i; j++){
                for (uint d = 0; d < 3; d++) dr[d] = Pos[j][d] - Pos[i][d];
            #ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #endif
                // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
                E_Potential[i] -= Mass[i] * Mass[j] / r;
            }
            for (uint j = i + 1; j < Size; j++){
                // Newton's third law: F_ij = -F_ji
                for (uint d = 0; d < 3; d++) dr[d] = -(Pos[i][d] - Pos[j][d]);
            #ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            #endif
                // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
                E_Potential[i] -= Mass[i] * Mass[j] / r;
            }
        }
#ifdef OPEN_MP
    }
#endif

    if (call_counter == 0)
        cout << __func__ << "...done!" << endl;
    call_counter++;
} // FUNCTION : Compute_PotentialEnergy_CPU

// Compute the total potential energy of the system
real ComputeTotalPotentialEnergy_CPU(real *E_Potential, const uint Size){
    static uint call_counter = 0;
    if (call_counter == 0) cout << "Which function is running?   " << __func__ << endl;
    real E_Potential_total = 0.0;
#ifdef SIMD
    // Use SIMD to speed up the summation
    #pragma omp simd reduction(+:E_Potential_total)
#endif
    // loop over all the particles
    for (uint i = 0; i < Size; i++){
        E_Potential_total += E_Potential[i];
    }
    cout << "Total Potential Energy: " << E_Potential_total << endl;

    if (call_counter == 0) cout << __func__ << "...done!" << endl;
    call_counter++;
    return E_Potential_total;
}

