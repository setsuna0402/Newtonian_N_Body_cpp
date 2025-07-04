/***********************************************************************
/
/  Compute potential energy of each particle by soften gravity (GPU version)
/  Note that Gravitational constant is always assumed to be 1 !
/  There are three GPU algorithms to compute the potential energy.
/  Compute_PotentialEnergy_GPU_SLOW: This is a slow version, 
/  which is not optimised for performance. 
/  (No per-thread register, no shared memory usage, etc.)
/  Compute_PotentialEnergy_GPU_FAST: This is faster than Compute_PotentialEnergy_GPU_SLOW, 
/  using per-thread register.
/  Compute_PotentialEnergy_GPU_SHARED: This is a shared memory version,
/  which use per-thread register and shared memory to optimise the performance.
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


// GPU function : Compute the potential for each particles
__global__ void Compute_PotentialEnergy_GPU_SLOW(real (*Pos)[3], real *Mass,
                                                 real *E_Potential, const uint Size){

    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0.0f, 0.0f, 0.0f};
    real r = 0.0f;
    real r3 = 0.0f;
    // get the target particle index of each thread
    const int i = blockDim.x*blockIdx.x + threadIdx.x;
    // calculate the potential energy for the i-th particle from all j particle
    // dr_ij : Position vector points to j-th particle from the i-th particle
    if (i < Size) {
        // initialise PE to zero
        E_Potential[i] = 0.0f;
        // loop over all the particles
        // Note E_Potential_i_i has to be a zero because there is no self gravity
        // Seperate to two loop to ensure the self gravity is not counted
        for (uint j=0; j<i; j++) {
            for (uint d=0; d<3; d++) dr[d] = Pos[j][d] - Pos[i][d];
#ifdef FLOAT8
            r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
#else
            r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
#endif
            r3 = r * r * r;
            // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
            E_Potential[i] -= Mass[i] * Mass[j] / r;
        }
        for (uint j=i+1; j<Size; j++) {
            // for (uint d=0; d<3; d++) dr[d] = -(Pos[i][d] - Pos[j][d]);
            for (uint d=0; d<3; d++) dr[d] = Pos[i][d] - Pos[j][d];
            // for (uint d=0; d<3; d++) dr[d] = Pos[j][d] - Pos[i][d];

#ifdef FLOAT8
            r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            // r3 = pow(r, 3.0);
#else
            r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            // r3 = powf(r, 3.0f);
#endif
            r3 = r * r * r;
            // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
            E_Potential[i] -= Mass[i] * Mass[j] / r;
        }
    }
    // printf("id = %i : E_Potential = %e \n", i, E_Potential[i]);
} // FUNCTION : Compute_PotentialEnergy_GPU_SLOW

// GPU function : Compute the PE (using per-thread register) for each particles
__global__ void Compute_PotentialEnergy_GPU_FAST(real (*Pos)[3], real *Mass,
                                                 real *E_Potential, const uint Size){
    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0.0f, 0.0f, 0.0f};
    real r = 0.0f;
    real r3 = 0.0f;
    // per-thread register to store the mass, PE and position of the i-th particle
    real mass_thread = 0.0f;
    real pos_thread[3] = {0.0f, 0.0f, 0.0f};
    real e_potential_thread = 0.0f;
    // get the target particle index of each thread
    const int i = blockDim.x*blockIdx.x + threadIdx.x;
    // calculate the potential energy for the i-th particle from all j particle
    // dr_ij : Position vector points to j-th particle from the i-th particle
    if (i < Size) {
        // load data from global memory to the **per-thread** registers (much faster)
        for (int d=0; d<3; d++) pos_thread[d] = Pos[i][d];
        mass_thread = Mass[i];  

        // loop over all the particles
        // Note E_Potential_i_i has to be a zero because there is no self gravity
        // Seperate to two loop to ensure the self gravity is not counted
        for (uint j=0; j<i; j++) {
            for (uint d=0; d<3; d++) dr[d] = Pos[j][d] - pos_thread[d];
#ifdef FLOAT8
            r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            // r3 = pow(r, 3.0f);
#else
            r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            // r3 = powf(r, 3.0f);
#endif
            r3 = r * r * r;
            // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
            e_potential_thread -= mass_thread * Mass[j] / r;
        }
        for (uint j=i+1; j<Size; j++) {
            // for (uint d=0; d<3; d++) dr[d] = -(pos_thread[d] - Pos[j][d]);
            for (uint d=0; d<3; d++) dr[d] = pos_thread[d] - Pos[j][d];
            // for (uint d=0; d<3; d++) dr[d] = Pos[j][d] - pos_thread[d];

#ifdef FLOAT8
            r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            // r3 = pow(r, 3.0);
#else
            r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
            // r3 = powf(r, 3.0f);
#endif
            r3 = r * r * r;
            // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
            e_potential_thread -= mass_thread * Mass[j] / r;
        }
        // load the PE from the per-thread register to the global memory
        E_Potential[i] = e_potential_thread;
    }
    // printf("id = %i : E_Potential = %e \n", i, E_Potential[i]);
} // FUNCTION : Compute_PotentialEnergy_GPU_FAST

// GPU function : Compute the PE (using shared memory) for each particles
__global__ void Compute_PotentialEnergy_GPU_SHARED(real (*Pos)[3], real *Mass,
                                               real *E_Potential, const uint Size){
    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0.0f, 0.0f, 0.0f};
    real r = 0.0f;
    real r3 = 0.0f;
    // per-thread register to store the mass, PE and position of the i-th particle
    real mass_thread = 0.0f;
    real pos_thread[3] = {0.0f, 0.0f, 0.0f};
    real e_potential_thread = 0.0f;
    uint j;
    int Threshold_Size = 0;
    int Threshold_sysmmety = 0;
    // get the target particle index of each thread
    const int i   = blockDim.x*blockIdx.x + threadIdx.x;
    const int tid = threadIdx.x;
    // declared shared-memory arrays
    __shared__ real share_Pos_x[GPU_BLOCK_SIZE];
    __shared__ real share_Pos_y[GPU_BLOCK_SIZE];
    __shared__ real share_Pos_z[GPU_BLOCK_SIZE];
    __shared__ real share_Mass [GPU_BLOCK_SIZE];
    // calculate the potential energy for the i-th particle from all j particle
    // dr_ij : Position vector points to j-th particle from the i-th particle
    // load data from global memory to the **per-thread** registers (much faster)
    if (i < Size) {
        for (int d=0; d<3; d++) pos_thread[d] = Pos[i][d];
        mass_thread = Mass[i];  
    }    
    for (int J_Base=0; J_Base<Size; J_Base+=GPU_BLOCK_SIZE) {
        // synchronise all threads before loading data
        __syncthreads();
        j = J_Base + tid;
        // load data from global memory to shared-memory
        if (j < Size) {
            share_Pos_x[tid] = Pos[j][0];
            share_Pos_y[tid] = Pos[j][1];
            share_Pos_z[tid] = Pos[j][2];
            share_Mass [tid] = Mass[j]  ;
        }
        __syncthreads(); // synchronise all threads to ensure all data have been loaded
        if (i >= Size) continue;  // No particles if i >= Size
        // Be careful the actual number of particles in the shared memory
        // Don't calculate non-exist particles
        // Below calcuate the actual number of particles in the shared memory
        Threshold_Size = min(GPU_BLOCK_SIZE, (Size - J_Base));
        if ((J_Base + Threshold_Size) <= i) {
            for (int k=0; k<Threshold_Size; k++) {
                // if (J_Base + Threshold_Size) == i, J_Base + k_max == i - 1
                // Ensure PE_i_i is zero!
                dr[0] = share_Pos_x[k] - pos_thread[0];
                dr[1] = share_Pos_y[k] - pos_thread[1];
                dr[2] = share_Pos_z[k] - pos_thread[2];
#ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = pow(r, 3.0f);
#else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = powf(r, 3.0f);
#endif
                r3 = r * r * r;
                // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
                e_potential_thread -= mass_thread * share_Mass[k] / r;
            }
        } else if (J_Base > i) {
            for (int k=0; k<Threshold_Size; k++) {
                dr[0] = pos_thread[0] - share_Pos_x[k];
                dr[1] = pos_thread[1] - share_Pos_y[k];
                dr[2] = pos_thread[2] - share_Pos_z[k];
#ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = pow(r, 3.0f);
#else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = powf(r, 3.0f);
#endif
                r3 = r * r * r;
                // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
                e_potential_thread -= mass_thread * share_Mass[k] / r;
            }
        } else {  // J_Base <= i < J_Base + Threshold_Size
            Threshold_sysmmety = i - J_Base;
            for (int k=0; k<Threshold_sysmmety; k++) {
                // When i == J_Base ==> Threshold_sysmmety == 0,
                // So this skips the i-th loop to avoid self-gravity
                dr[0] = share_Pos_x[k] - pos_thread[0];
                dr[1] = share_Pos_y[k] - pos_thread[1];
                dr[2] = share_Pos_z[k] - pos_thread[2];
#ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = pow(r, 3.0f);
#else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = powf(r, 3.0f);
#endif
                r3 = r * r * r;
                // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
                e_potential_thread -= mass_thread * share_Mass[k] / r;
            }
            for (int k=Threshold_sysmmety+1; k<Threshold_Size; k++) {
                dr[0] = pos_thread[0] - share_Pos_x[k];
                dr[1] = pos_thread[1] - share_Pos_y[k];
                dr[2] = pos_thread[2] - share_Pos_z[k];
#ifdef FLOAT8
                r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = pow(r, 3.0f);
#else
                r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                // r3 = powf(r, 3.0f);
#endif
                r3 = r * r * r;
                // accumulate the PE (PE_ij = Mass_i * Mass_j / |dr_ij|)
                e_potential_thread -= mass_thread * share_Mass[k] / r;
            }
        }
    } // for (int J_Base=0; J_Base<Size; J_Base+=GPU_BLOCK_SIZE)
    // load the PE from the per-thread register to the global memory
    if (i < Size) E_Potential[i] = e_potential_thread;

    // printf("id = %i : E_Potential = %e \n", i, E_Potential[i]);
} // FUNCTION : Compute_PotentialEnergy_GPU_SHARED

