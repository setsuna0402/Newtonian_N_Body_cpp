/***********************************************************************
/
/  Soften gravity N Body simulation (GPU): 
/  evolve the direct N-Body system by GPU
/
/  written by: KH 
/  date:       2025/6/19
/  modified1:
/
/  Note: 
/   Unit system:
/   Box size = 1
/   Gravitational Constant assumes to 1. 
/   So mass unit is determinated by distance and time units.
/   User need to make sure the unit conversion factors are correct:
/   G_Length_CodeUnits: convert distance from code unit to meter
/   G_Time_CodeUnits  : convert time from code unit to second
/   G_Mass_CodeUnits  : convert mass from code unit to kg
/   G_Velocity_CodeUnits: convert velocity from code unit to meter/second
/   G_Energy_CodeUnits: convert energy from code unit to Joule
/   Coordinate order: [x, y, z]
/   Use Leapfrog (kick-drift-kick) to update the position and velocity
/   Use the soften gravity method to calculate the acceleration
/   C++ - GPU parallelisation (CUDA)
/   There is only one GPU stream in this code, so GPU functions are called sequentially
/
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"
using namespace std;
// Declare global constant
// Gravitational constant in mks unit m^3/(kg⋅s^2)
const real Gc_GravConst_mks = 6.6743E-11       ;  // Gravitational constant in mks unit
const real Gc_GravConst_cgs = 6.6743E-8        ;  // cm^3/(g⋅s^2)
const real Gc_sqrt_pi       = 1.77245385091    ;  // sqrt(pi)
const real Gc_Mpc_to_km     = 3.08567758128E+19;  // 1Mpc = 3.08567758128E+19 km
const real Gc_Mpc_to_m      = 3.08567758128E+22;  // 1Mpc = 3.08567758128E+22 m
const real Gc_Mpc_to_cm     = 3.08567758128E+24;  // 1Mpc = 3.08567758128E+24 cm
const real Gc_Boxsize       = BOXSIZE_CODE     ;  // in code unit
const real Gc_End_time      = End_Time         ;  // Endtime of simulation (code unit)

// Declare global variables for unit conversion
real G_Length_CodeUnits   ;  // Convert distance from code unit to meter
real G_Time_CodeUnits     ;  // Convert time from code unit to second
real G_Mass_CodeUnits     ;  // Convert mass from code unit to kg
real G_Velocity_CodeUnits ;  // Convert velocity from code unit to meter/second
real G_Energy_CodeUnits   ;  // Convert energy from code unit to Joule

int main( int argc, char **argv ){
    
    const uint RSeed = 1234; // random seed
    srand(RSeed);
    clock_t start_time, end_time;

    // Variabls for calling GPU functions
    uint N_GPU_BLOCK = (N_PARTICLE + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE; // Number of blocks

    // Declare variables
    real dt_host      = 0.0; // Time step in code unit
    real dt_next_host = 0.0;
    uint num_round    = 0  ; // Number of round for evolution
    real time_host    = 0.0; // Total simulation time in code unit
    uint N_step_total = 0  ; // Total number of steps
    real PE_total     = 0.0; // Total potential energy of the system
    real KE_total     = 0.0; // Total kinetic energy of the system
    real Energy_total = 0.0; // Total energy of the system
#ifdef FLOAT8
    const real const_half = 0.5; // Half factor for velocity update
#else
    const real const_half = 0.5f; // Half factor for velocity update
#endif

    // Declare pointers (for arrays) in host (CPU)
    real  *h_Mass        = NULL;  // Mass of each particle
    real  *h_E_Potential = NULL;  // Potential energy of each particle
    real  *h_E_Kinetic   = NULL;  // Kinetic energy of each particle
    real (*h_Pos)[3]     = NULL;  // Position of each particle
    real (*h_Vel)[3]     = NULL;  // Velocity of each particle
    real (*h_Acc)[3]     = NULL;  // Acceleration of each particle
    
    // Declare pointers (for arrays) in device (GPU)
    real  *d_Mass        = NULL;  // Mass of each particle
    real  *d_E_Potential = NULL;  // Potential energy of each particle
    real  *d_E_Kinetic   = NULL;  // Kinetic energy of each particle
    real (*d_Pos)[3]     = NULL;  // Position of each particle
    real (*d_Vel)[3]     = NULL;  // Velocity of each particle
    real (*d_Acc)[3]     = NULL;  // Acceleration of each particle
    real  *d_sum_part_pe = NULL;  // Use to store the local sum of each GPU block
    real  *d_sum_all_pe  = NULL;  // Use to store the total sum of all GPU blocks
    real  *d_sum_part_ke = NULL;  // Use to store the local sum of each GPU block
    real  *d_sum_all_ke  = NULL;  // Use to store the total sum of all GPU blocks

    // Allocate memory in host (CPU)
    h_Mass        = new real[N_PARTICLE]   ;
    h_E_Potential = new real[N_PARTICLE]   ;
    h_E_Kinetic   = new real[N_PARTICLE]   ;
    h_Pos         = new real[N_PARTICLE][3];
    h_Vel         = new real[N_PARTICLE][3];
    h_Acc         = new real[N_PARTICLE][3];

    // Allocate memory in device (GPU)
    // Relevant to the particle number
    CHECK_CUDA(cudaMalloc(&d_Mass       ,     N_PARTICLE * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_E_Potential,     N_PARTICLE * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_E_Kinetic  ,     N_PARTICLE * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_Pos        , 3 * N_PARTICLE * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_Vel        , 3 * N_PARTICLE * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_Acc        , 3 * N_PARTICLE * sizeof(real)));
    // Relevant to the GPU block number (Use to do reduction sum)
    CHECK_CUDA(cudaMalloc(&d_sum_part_pe, N_GPU_BLOCK * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_sum_part_ke, N_GPU_BLOCK * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_sum_all_pe , 1 * sizeof(real)));
    CHECK_CUDA(cudaMalloc(&d_sum_all_ke , 1 * sizeof(real)));

/*
    if (MetaParameterReader() != EXIT_SUCCESS) {
        fprintf(stderr, "Func MetaParameterReader() has some issues,\n");
        fprintf(stderr, "Please check the meta parameter file!\n");
        fprintf(stderr, "Terminating!\n");
        return EXIT_FAILURE;
    }
*/
    // Set all code units to one
    if (MetaParameterUnity() != 0) {
        fprintf(stderr, "Func MetaParameterUnity() has some issues.\n");
        fprintf(stderr, "Terminating!\n");
        return EXIT_FAILURE;
    }
    // Transfer the global varsiabls to device 
    if (PassGlobleVariabletoDevice() != EXIT_SUCCESS) {
        fprintf(stdout, "PassGlobleVariabletoDevice() failed, terminate the programe.\n");
        return EXIT_FAILURE;
    }

    start_time = clock();
    CreateRandomIC(0.25 * BOXSIZE_CODE, h_Pos, h_Vel, h_Mass, N_PARTICLE);
    end_time = clock();
    printf("CreateRandomIC() time cost: %f seconds.\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Copy the initial data from host to device
    CHECK_CUDA(cudaMemcpy(d_Mass , h_Mass,     N_PARTICLE * sizeof(real), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_Pos  , h_Pos , 3 * N_PARTICLE * sizeof(real), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_Vel  , h_Vel , 3 * N_PARTICLE * sizeof(real), cudaMemcpyHostToDevice));

    // Compute the initial acceleration

    start_time = clock();
#if defined(GPU_SLOW)
    Compute_ACC_GPU_SLOW <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
#elif defined(GPU_FAST)
    Compute_ACC_GPU_FAST <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
#elif defined(GPU_SHARED)
    Compute_ACC_GPU_SHARED <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
#else
    // Error: No GPU acceleration method defined
    fprintf(stderr, "Error: No GPU acceleration method defined.\n");
    fprintf(stderr, "Please define GPU_SLOW, GPU_FAST or GPU_SHARED in macros_and_parameters.h.\n");
    return EXIT_FAILURE;
#endif 
    // Check for errors in the kernel launch
    CHECK_CUDA(cudaGetLastError()     );
    // Copy the acceleration data from device to host
    CHECK_CUDA(cudaMemcpy(h_Acc, d_Acc, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));

    end_time = clock();
    printf("How long did we take for computing Acc? %f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Compute the initial time step
    Safest_TimeStep_CPU(dt_host);
    if (dt_host >= End_Time) dt_host = 0.1 * End_Time;
    // printf("time_host = %2.6f, dt = %2.6e\n", time_host, dt_host);
    N_step_total = (uint)(End_Time / dt_host) + 1; // Total number of steps
    if (N_step_total > N_STEP_MAX) {
        fprintf(stderr, "N_step_total = %d exceeds N_STEP_MAX = %d.\n", 
                N_step_total, N_STEP_MAX);
        fprintf(stderr, "Please reduce End_Time in macros_and_parameters.h.\n");
        return EXIT_FAILURE;
    }
    printf("End time = %f, N_step_total = %d, dt = %2.6e\n",
            End_Time, N_step_total, dt_host);
    
    start_time = clock();
    // Compute the initial potential energy
#if defined(GPU_SLOW)
    Compute_PotentialEnergy_GPU_SLOW <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
#elif defined(GPU_FAST)
    Compute_PotentialEnergy_GPU_FAST <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
#elif defined(GPU_SHARED)
    Compute_PotentialEnergy_GPU_SHARED <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
#endif
    CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
    // Initialise the arrays to zero
    CHECK_CUDA(cudaMemset(d_sum_part_pe, 0, N_GPU_BLOCK)); // Set all elements to 0
    CHECK_CUDA(cudaMemset(d_sum_all_pe , 0, 1          )); // Set all elements to 0
    // Perform reduction sum for potential energy
    ReductionSum_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_E_Potential, d_sum_part_pe, N_PARTICLE);
    ReductionSum_GPU <<< 1, GPU_BLOCK_SIZE >>> (d_sum_part_pe, d_sum_all_pe, N_GPU_BLOCK);
    CHECK_CUDA(cudaGetLastError());
    // Copy the result back to host
    CHECK_CUDA(cudaMemcpy(&PE_total, d_sum_all_pe, sizeof(real), cudaMemcpyDeviceToHost));
    // the potential energy is divided by 2, because we court all the pairs twice
#ifdef FLOAT8
    PE_total = 0.5 * PE_total;
#else
    PE_total = 0.5f * PE_total;
#endif
    end_time = clock();
    printf("How long did we take for computing PE? %f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Compute the initial kinetic energy
    start_time = clock();
    Compute_KineticEnergy_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Vel, d_Mass, d_E_Kinetic, N_PARTICLE);
    CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
    // Initialise the arrays to zero
    CHECK_CUDA(cudaMemset(d_sum_part_ke, 0, N_GPU_BLOCK)); // Set all elements to 0
    CHECK_CUDA(cudaMemset(d_sum_all_ke , 0, 1          )); // Set all elements to 0
    // Perform reduction sum for kinetic energy
    ReductionSum_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_E_Kinetic, d_sum_part_ke, N_PARTICLE);
    ReductionSum_GPU <<< 1, GPU_BLOCK_SIZE >>> (d_sum_part_ke, d_sum_all_ke, N_GPU_BLOCK);
    CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
    // Copy the result back to host
    CHECK_CUDA(cudaMemcpy(&KE_total, d_sum_all_ke, sizeof(real), cudaMemcpyDeviceToHost));
    end_time = clock();
    printf("How long did we take for computing KE? %f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);
    
    // Print the initial total energy
    Energy_total = PE_total + KE_total;
    printf("Initial Total Energy = %2.6e\n", Energy_total);
    // Copy the data from device to host
    CHECK_CUDA(cudaMemcpy(h_E_Potential, d_E_Potential, N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(h_E_Kinetic  , d_E_Kinetic  , N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    
    // check sum of potential and kinetic energy
    real PE_sum = 0.0, KE_sum = 0.0;
    for (uint i = 0; i < N_PARTICLE; i++) {
        PE_sum += h_E_Potential[i];
        KE_sum += h_E_Kinetic[i];
    }
    PE_sum = 0.5f * PE_sum; // Divide by 2 because we count all pairs twice
    printf("Sum (GPU) of potential energy = %2.6e, sum of kinetic energy = %2.6e\n",
            PE_total, KE_total);
    printf("Sum (CPU) of potential energy = %2.6e, sum of kinetic energy = %2.6e\n",
           PE_sum, KE_sum);
    // Save the initial snapshot
    SaveDataHDF5(h_Pos, h_Vel, h_Acc, h_Mass, h_E_Potential, h_E_Kinetic,
                 N_PARTICLE, num_round, dt_host, time_host);

    time_host = 0.0f; // Make sure time_host starts from 0.0
    cudaDeviceSynchronize();    // Synchronize the device before starting the evolution loop
    start_time = clock();
    // Start the evolution loop
    for (uint round = 1; round <= N_step_total; round++)
    {
        num_round = round;
        // Compute acc(t_i) using Pos(t_i)
    #if defined(GPU_SLOW)
        Compute_ACC_GPU_SLOW <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
    #elif defined(GPU_FAST)
        Compute_ACC_GPU_FAST <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
    #elif defined(GPU_SHARED)
        Compute_ACC_GPU_SHARED <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
    #endif
        CHECK_CUDA(cudaGetLastError());  
        // Update Velocity to t_i+1/2 (Kick)
        UpdateVelocity_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (dt_host * const_half, d_Vel, d_Acc, N_PARTICLE);
        CHECK_CUDA(cudaGetLastError());
        // Update Position to t_i+1 (Drift)
        UpdatePosition_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (dt_host, d_Pos, d_Vel, N_PARTICLE);
        CHECK_CUDA(cudaGetLastError());
        // Compute acc(t_i+1) using Pos(t_i+1)
    #if defined(GPU_SLOW)
        Compute_ACC_GPU_SLOW <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
    #elif defined(GPU_FAST)
        Compute_ACC_GPU_FAST <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
    #elif defined(GPU_SHARED)
        Compute_ACC_GPU_SHARED <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Acc, d_Mass, N_PARTICLE);
    #endif
        CHECK_CUDA(cudaGetLastError());
        // Update Velocity to t_i+1 (Kick)
        UpdateVelocity_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (dt_host * const_half, d_Vel, d_Acc, N_PARTICLE);
        CHECK_CUDA(cudaGetLastError());
        // Update simulation time
        time_host += dt_host;
        // check if the time has reached the end time
        if (time_host >= End_Time){
            printf("time_host = %2.6f, dt = %2.6e\n", time_host, dt_host);
            printf("Code time has reached the End_Time, stop now.\n");
            cudaDeviceSynchronize();    // Ensure all GPU tasks are completed
            break;
        }
        // check if we save a snapshot in this round
        if (round % N_DUMP_STEP == 0){
            cudaDeviceSynchronize();    // Synchronize the device
            printf("Total round = %d, current round = %d\n", N_step_total, round);
            printf("Saving a snapshot at time_host = %2.6f, dt = %2.6e\n",
                    time_host, dt_host);
            // Compute the potential energy
        #if defined(GPU_SLOW)
            Compute_PotentialEnergy_GPU_SLOW <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
        #elif defined(GPU_FAST)
            Compute_PotentialEnergy_GPU_FAST <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
        #elif defined(GPU_SHARED)
            Compute_PotentialEnergy_GPU_SHARED <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
        #endif
            CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
            // Initialise the arrays to zero
            CHECK_CUDA(cudaMemset(d_sum_part_pe, 0, N_GPU_BLOCK)); // Set all elements to 0
            CHECK_CUDA(cudaMemset(d_sum_all_pe , 0, 1          )); // Set all elements to 0
            // Perform reduction sum for potential energy
            ReductionSum_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_E_Potential, d_sum_part_pe, N_PARTICLE);
            ReductionSum_GPU <<< 1, GPU_BLOCK_SIZE >>> (d_sum_part_pe, d_sum_all_pe, N_GPU_BLOCK);
            CHECK_CUDA(cudaGetLastError());
            // Copy the result back to host
            CHECK_CUDA(cudaMemcpy(&PE_total, d_sum_all_pe, sizeof(real), cudaMemcpyDeviceToHost));

            // Compute the kinetic energy
            Compute_KineticEnergy_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Vel, d_Mass, d_E_Kinetic, N_PARTICLE);
            CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
            // Initialise the arrays to zero
            CHECK_CUDA(cudaMemset(d_sum_part_ke, 0, N_GPU_BLOCK)); // Set all elements to 0
            CHECK_CUDA(cudaMemset(d_sum_all_ke , 0, 1          )); // Set all elements to 0
            // Perform reduction sum for kinetic energy
            ReductionSum_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_E_Kinetic, d_sum_part_ke, N_PARTICLE);
            ReductionSum_GPU <<< 1, GPU_BLOCK_SIZE >>> (d_sum_part_ke, d_sum_all_ke, N_GPU_BLOCK);
            CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
            // Copy the result back to host
            CHECK_CUDA(cudaMemcpy(&KE_total, d_sum_all_ke, sizeof(real), cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();        // Synchronize the device
            // the potential energy is divided by 2, because we court all the pairs twice
        #ifdef FLOAT8
            PE_total = 0.5 * PE_total;
        #else
            PE_total = 0.5f * PE_total;
        #endif
            // Compute the total energy
            Energy_total = PE_total + KE_total;
            // Print the total energy
            printf("Total Energy = %2.6e\n", Energy_total);
            // Copy the data from device to host
            CHECK_CUDA(cudaMemcpy(h_Pos, d_Pos, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
            CHECK_CUDA(cudaMemcpy(h_Vel, d_Vel, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
            CHECK_CUDA(cudaMemcpy(h_Acc, d_Acc, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
            CHECK_CUDA(cudaMemcpy(h_E_Potential, d_E_Potential, N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
            CHECK_CUDA(cudaMemcpy(h_E_Kinetic  , d_E_Kinetic  , N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
            // Save the snapshot
            SaveDataHDF5(h_Pos, h_Vel, h_Acc, h_Mass, h_E_Potential, h_E_Kinetic,
                         N_PARTICLE, num_round, dt_host, time_host);
            printf("Snapshot saved successfully!\n");
        }
    }
    cudaDeviceSynchronize();    // Ensure the evolution loop is completed
    // Compute the potential energy
#if defined(GPU_SLOW)
    Compute_PotentialEnergy_GPU_SLOW <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
#elif defined(GPU_FAST)
    Compute_PotentialEnergy_GPU_FAST <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
#elif defined(GPU_SHARED)
    Compute_PotentialEnergy_GPU_SHARED <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Pos, d_Mass, d_E_Potential, N_PARTICLE);
#endif
    CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
    // Initialise the arrays to zero
    CHECK_CUDA(cudaMemset(d_sum_part_pe, 0, N_GPU_BLOCK)); // Set all elements to 0
    CHECK_CUDA(cudaMemset(d_sum_all_pe , 0, 1          )); // Set all elements to 0
    // Perform reduction sum for potential energy
    ReductionSum_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_E_Potential, d_sum_part_pe, N_PARTICLE);
    ReductionSum_GPU <<< 1, GPU_BLOCK_SIZE >>> (d_sum_part_pe, d_sum_all_pe, N_GPU_BLOCK);
    CHECK_CUDA(cudaGetLastError());
    // Copy the result back to host
    CHECK_CUDA(cudaMemcpy(&PE_total, d_sum_all_pe, sizeof(real), cudaMemcpyDeviceToHost));
    
    // Compute the kinetic energy
    Compute_KineticEnergy_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_Vel, d_Mass, d_E_Kinetic, N_PARTICLE);
    CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
    // Initialise the arrays to zero
    CHECK_CUDA(cudaMemset(d_sum_part_ke, 0, N_GPU_BLOCK)); // Set all elements to 0
    CHECK_CUDA(cudaMemset(d_sum_all_ke , 0, 1          )); // Set all elements to 0
    // Perform reduction sum for kinetic energy
    ReductionSum_GPU <<< N_GPU_BLOCK, GPU_BLOCK_SIZE >>> (d_E_Kinetic, d_sum_part_ke, N_PARTICLE);
    ReductionSum_GPU <<< 1, GPU_BLOCK_SIZE >>> (d_sum_part_ke, d_sum_all_ke, N_GPU_BLOCK);
    CHECK_CUDA(cudaGetLastError()); // Check for errors in the kernel launch
    // Copy the result back to host
    CHECK_CUDA(cudaMemcpy(&KE_total, d_sum_all_ke, sizeof(real), cudaMemcpyDeviceToHost));
            
    cudaDeviceSynchronize(); // Ensure we get PE_total and KE_total before printing
    // the potential energy is divided by 2, because we court all the pairs twice
#ifdef FLOAT8
    PE_total = 0.5 * PE_total;
#else
    PE_total = 0.5f * PE_total;
#endif
    // Compute the total energy
    Energy_total = PE_total + KE_total;
    // Print the total energy
    printf("Final Total Energy = %2.6e at t = %2.6e\n", Energy_total, time_host);

    // Copy the data from device to host
    CHECK_CUDA(cudaMemcpy(h_Pos, d_Pos, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(h_Vel, d_Vel, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(h_Acc, d_Acc, 3 * N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(h_E_Potential, d_E_Potential, N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(h_E_Kinetic  , d_E_Kinetic  , N_PARTICLE * sizeof(real), cudaMemcpyDeviceToHost));
    SaveDataHDF5(h_Pos, h_Vel, h_Acc, h_Mass, h_E_Potential, h_E_Kinetic,
        N_PARTICLE, num_round, dt_host, time_host);
    printf("The final snapshot saved successfully!\n");

    end_time = clock();
    printf("This simulation costs %f seconds with GPU!\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Release memory in device (GPU)
    CHECK_CUDA(cudaFree(d_Mass       ));
    CHECK_CUDA(cudaFree(d_E_Potential));
    CHECK_CUDA(cudaFree(d_E_Kinetic  ));
    CHECK_CUDA(cudaFree(d_Pos        ));
    CHECK_CUDA(cudaFree(d_Vel        ));
    CHECK_CUDA(cudaFree(d_Acc        ));
    CHECK_CUDA(cudaFree(d_sum_part_pe));
    CHECK_CUDA(cudaFree(d_sum_part_ke));
    CHECK_CUDA(cudaFree(d_sum_all_pe ));
    CHECK_CUDA(cudaFree(d_sum_all_ke ));

    // Release memory in host (CPU)
    delete[] h_Mass;
    delete[] h_E_Potential;
    delete[] h_E_Kinetic;
    delete[] h_Pos;
    delete[] h_Vel;
    delete[] h_Acc;
    return EXIT_SUCCESS;        
}   // End of main function
