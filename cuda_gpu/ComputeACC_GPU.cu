/***********************************************************************
/
/  Compute acceleration of each particle by soften gravity (GPU version)
/  Note that Gravitational constant is always assumed to be 1 !
/  There are three GPU algorithms to compute the acceleration.
/  Compute_ACC_GPU_SLOW: This is a slow version, which is not optimised for performance.
/  (No per-thread register, no shared memory usage, etc.)
/  Compute_ACC_GPU_FAST: This is a faster than Compute_ACC_GPU_SLOW,
/  using per-thread register.
/  Compute_ACC_GPU_SHARED: This is a shared memory version,
/  which use per-thread register and shared memory to optimise the performance.
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

// GPU function : Compute the acceleration for each particles by soften gravity
__global__ void Compute_ACC_GPU_SLOW(real (*Pos)[3], real (*Acc)[3],
                                     real *Mass, const uint Size){

    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0.0f, 0.0f, 0.0f};
    real r = 0.0f;
    real r3 = 0.0f;
    // get the target particle index of each thread
    const int i = blockDim.x*blockIdx.x + threadIdx.x;
    // calculate the acceleration for the i-th particle from all j particle
    //  F_ij : Force on the i-th particle caused by the j particle
    if (i < Size) {
        // initialise acc to zero
        // for (uint d=0; d<3; d++) Acc[i][d] = 0.0f; 
        Acc[i][0] = 0.0f;
        Acc[i][1] = 0.0f;
        Acc[i][2] = 0.0f;
        // loop over all the particles
        // Note Acc_i_i has to be a zero because there is no self gravity
        // Seperate to two loop to follow Newton's second law
        for (uint j=0; j<i; j++) {
            for (uint d=0; d<3; d++) dr[d] = Pos[j][d] - Pos[i][d];
#ifdef FLOAT8
            r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
#else
            r = sqrtf(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
#endif
            r3 = r * r * r;
            // accumulate the acceleration (acc = Mass_j r/|r|^3)
            for (uint d=0; d<3; d++) Acc[i][d] += Mass[j] * dr[d] / r3;
        }
        for (uint j=i+1; j<Size; j++){
            // I think the first line can eusure the newton's third law (numerically)
            // But it seems like compiller optimised this line
            // and it becomes the same as the third line
            // so, base on the KISS principle, I use the second line in all GPU functions
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
            // accumulate the acceleration (acc = Mass_j r/|r|^3)
            for (uint d=0; d<3; d++) Acc[i][d] -= Mass[j] * dr[d] / r3;
        }
    }
    // printf("id = %i : Acc[0] = %e, Acc[1] = %e, Acc[2] = %e \n", i, Acc[i][0], Acc[i][1], Acc[i][2]);
} // FUNCTION : Compute_ACC_GPU_SLOW

// GPU function : Compute the acceleration (using per-thread register) for each particles
__global__ void Compute_ACC_GPU_FAST(real (*Pos)[3], real (*Acc)[3],
                                     real *Mass, const uint Size){
    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0.0f, 0.0f, 0.0f};
    real r = 0.0f;
    real r3 = 0.0f;
    // per-thread register to store the mass, acceleration and position of the i-th particle
    real mass_thread = 0.0f;
    real acc_thread[3] = {0.0f, 0.0f, 0.0f};
    real pos_thread[3] = {0.0f, 0.0f, 0.0f};
    // get the target particle index of each thread
    const int i = blockDim.x*blockIdx.x + threadIdx.x;
    // calculate the acceleration for the i-th particle from all j particle
    //  F_ij : Force on the i-th particle caused by the j particle
    if (i < Size) {
        // load data from global memory to the **per-thread** registers (much faster)
        for (int d=0; d<3; d++)    pos_thread[d] = Pos[i][d];
        // mass_thread = Mass[i];  // This is not used for computing acc

        // loop over all the particles
        // Note Acc_i_i has to be a zero because there is no self gravity
        // Seperate to two loop to follow Newton's second law
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
            // accumulate the acceleration (acc = Mass_j r/|r|^3)
            for (uint d=0; d<3; d++) acc_thread[d] += Mass[j] * dr[d] / r3;
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
            // accumulate the acceleration (acc = Mass_j r/|r|^3)
            for (uint d=0; d<3; d++) acc_thread[d] -= Mass[j] * dr[d] / r3;
        }
        // load the acceleration from the per-thread register to the global memory
        for (uint d=0; d<3; d++) Acc[i][d] = acc_thread[d];
    }
    // printf("id = %i : Acc[0] = %e, Acc[1] = %e, Acc[2] = %e \n", i, Acc[i][0], Acc[i][1], Acc[i][2]);
} // FUNCTION : Compute_ACC_GPU_FAST

// GPU function : Compute the acceleration (using shared memory) for each particles
__global__ void Compute_ACC_GPU_SHARED(real (*Pos)[3], real (*Acc)[3],
                                       real *Mass, const uint Size){
    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    real dr[3] = {0.0f, 0.0f, 0.0f};
    real r = 0.0f;
    real r3 = 0.0f;
    // per-thread register to store the mass, acceleration and position of the i-th particle
    real mass_thread = 0.0f;
    real acc_thread[3] = {0.0f, 0.0f, 0.0f};
    real pos_thread[3] = {0.0f, 0.0f, 0.0f};
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
    // calculate the acceleration for the i-th particle from all j particle
    //  F_ij : Force on the i-th particle caused by the j particle
    // load data from global memory to the **per-thread** registers (much faster)
    if (i < Size) {
        for (int d=0; d<3; d++)    pos_thread[d] = Pos[i][d];
        // mass_thread = Mass[i];  // This is not used for computing acc
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
                // Ensure Acc_i_i is zero!
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
                // accumulate the acceleration (acc = Mass_j r/|r|^3)
                for (uint d=0; d<3; d++) acc_thread[d] += share_Mass[k] * dr[d] / r3;
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
                // accumulate the acceleration (acc = Mass_j r/|r|^3)
                for (uint d=0; d<3; d++) acc_thread[d] -= share_Mass[k] * dr[d] / r3;
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
                // accumulate the acceleration (acc = Mass_j r/|r|^3)
                for (uint d=0; d<3; d++) acc_thread[d] += share_Mass[k] * dr[d] / r3;
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
                // accumulate the acceleration (acc = Mass_j r/|r|^3)
                for (uint d=0; d<3; d++) acc_thread[d] -= share_Mass[k] * dr[d] / r3;
            }
        }
    } // for (int J_Base=0; J_Base<Size; J_Base+=GPU_BLOCK_SIZE)
    if (i < Size) for (uint d=0; d<3; d++) Acc[i][d] = acc_thread[d];
    // printf("id = %i : Acc[0] = %e, Acc[1] = %e, Acc[2] = %e \n", i, Acc[i][0], Acc[i][1], Acc[i][2]);
} // FUNCTION : Compute_ACC_GPU_SHARED
