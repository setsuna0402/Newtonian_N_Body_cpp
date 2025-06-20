/***********************************************************************
/
/  Compute kinetic energy of each particle and the total kinetic energy 
/    of the system
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Mass (Dynamical array) : Mass of each particle (Size: N)
/  E_Kinetic (Dynamical array) : Kinetic energy of each particle (Size: N)
/  Size (const uint ) : Number of particles


************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Compute the kinetic energy for each particles by soften gravity
real Compute_KineticEnergy_CPU(real (*Vel)[3], real *Mass, real *E_Kinetic,
                               const uint Size){
    static uint call_counter = 0;
    if (call_counter == 0) cout << "Which function is running?   " << __func__ << endl;
    real E_Kinetic_total = 0.0; // total kinetic energy of the system

    // Initialise the kinetic energy of each particle to zero
#ifdef FLOAT8
    for(uint i = 0; i < Size; i++) E_Kinetic[i] = 0.0;
#else
    for(uint i = 0; i < Size; i++) E_Kinetic[i] = 0.0f;
#endif
    // calculate the kinetic energy for the i-th particle
#ifdef OPEN_MP
    // Time cost of this loop many be small, so parallelisation may not be beneficial.
    // Use OpenMP to parallelise the for loop
    #pragma omp parallel
    {
        #pragma omp for
#endif
        for(uint i = 0; i < Size; i++){
            // KE_i = 1/2 * Mass_i * |Vel_i|^2
            // |Vel_i|^2 = Vel_i[0]^2 + Vel_i[1]^2 + Vel_i[2]^2
        #ifdef FLOAT8
            E_Kinetic[i] = 0.5 * Mass[i] * (Vel[i][0] * Vel[i][0] +
                                            Vel[i][1] * Vel[i][1] +
                                            Vel[i][2] * Vel[i][2]  );
        #else
            E_Kinetic[i] = 0.5f * Mass[i] * (Vel[i][0] * Vel[i][0] +
                                             Vel[i][1] * Vel[i][1] +
                                             Vel[i][2] * Vel[i][2]  );
        #endif
        }
#ifdef OPEN_MP
    }
#endif

    // calculate the total kinetic energy of the system
#ifdef SIMD
    // Use SIMD to parallelise the for loop
    #pragma omp simd reduction(+:E_Kinetic_total)
#endif
    for(uint i = 0; i < Size; i++) E_Kinetic_total += E_Kinetic[i];


    if (call_counter == 0)
        cout << __func__ << "...done!" << endl;
    call_counter++;
} // FUNCTION : Compute_KineticEnergy_CPU

// Compute the total kinetic energy of the system
real ComputeTotalKineticEnergy_CPU(real *E_Kinetic, const uint Size){
    static uint call_counter = 0;
    if (call_counter == 0) cout << "Which function is running?   " << __func__ << endl;
    real E_Kinetic_total = 0.0;
#ifdef SIMD
    // Use SIMD to speed up the summation
    #pragma omp simd reduction(+:E_Kinetic_total)
#endif
    // loop over all the particles
    for (uint i = 0; i < Size; i++){
        E_Kinetic_total += E_Kinetic[i];
    }
#ifdef DEBUG
    printf("Total Kinetic Energy: %f\n", E_Kinetic_total);
#endif
    if (call_counter == 0) cout << __func__ << "...done!" << endl;
    call_counter++;
    return E_Kinetic_total;
}// FUNCTION : ComputeTotalKineticEnergy_CPU
