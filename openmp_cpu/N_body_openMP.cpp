/***********************************************************************
/
/  SimpleNBody_CPU: use CPU to evolve the direct N-Body system
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
/   Use the softern gravity method to calculate the acceleration
/   C++ - multiple threads (OpenMP)
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
    // Initialise OpenMP
#ifdef OPEN_MP
    int nProcessors = omp_get_max_threads();
    if (nProcessors < OPENMP_NUM_THREAD) {
        fprintf(stderr, "Warning: nProcessors = %d is less than OPENMP_NUM_THREAD = %d.\n", 
                nProcessors, OPENMP_NUM_THREAD);
        fprintf(stderr, "Using nProcessors = %d instead.\n", nProcessors);
        omp_set_num_threads(nProcessors) ;    // Set how many thread we can use
    }
    else{
        omp_set_num_threads(OPENMP_NUM_THREAD);  // Set how many thread we can use
    }
    // omp_set_num_threads(nProcessors) ;    // Set how many thread we can use
    nProcessors = omp_get_max_threads();
    cout << "omp_get_max_threads = " << nProcessors << endl;
    double itime, ftime, exec_time;
#endif
    const uint RSeed = 1234; // random seed
    srand(RSeed);
    clock_t start_time, end_time;

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

    // Declare data arrays
    real  *h_Mass        = NULL;  // Mass of each particle
    real  *h_E_Potential = NULL;  // Potential energy of each particle
    real  *h_E_Kinetic   = NULL;  // Kinetic energy of each particle
    real (*h_Pos)[3]     = NULL;  // Position of each particle
    real (*h_Vel)[3]     = NULL;  // Velocity of each particle
    real (*h_Acc)[3]     = NULL;  // Acceleration of each particle
    // Allocate memory for particles
    h_Mass        = new real[N_PARTICLE]   ;
    h_E_Potential = new real[N_PARTICLE]   ;
    h_E_Kinetic   = new real[N_PARTICLE]   ;
    h_Pos         = new real[N_PARTICLE][3];
    h_Vel         = new real[N_PARTICLE][3];
    h_Acc         = new real[N_PARTICLE][3];

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

    start_time = clock();
    CreateRandomIC(0.25 * BOXSIZE_CODE, h_Pos, h_Vel, h_Mass, N_PARTICLE);
    end_time = clock();
    printf("CreateRandomIC() time cost: %f seconds.\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Compute the initial acceleration
#ifdef OPEN_MP
    itime = omp_get_wtime();
#else
    start_time = clock();
#endif
    Compute_ACC_CPU(h_Pos, h_Acc, h_Mass, N_PARTICLE);
#ifdef OPEN_MP
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("How long did we take for computing Acc? %f seconds\n", exec_time);
#else
    end_time = clock();
    printf("How long did we take for computing Acc? %f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif
    // Compute the initial time step
    Safest_TimeStep_CPU(dt_host);
    if (dt_host >= End_Time) dt_host = 0.1 * End_Time;
    // printf("time_host = %2.6f, dt = %2.6e\n", time_host, dt_host);
    N_step_total = (uint)(End_Time / dt_host); // Total number of steps
    if (N_step_total > N_STEP_MAX) {
        fprintf(stderr, "N_step_total = %d exceeds N_STEP_MAX = %d.\n", 
                N_step_total, N_STEP_MAX);
        fprintf(stderr, "Please reduce End_Time in macros_and_parameters.h.\n");
        return EXIT_FAILURE;
    }
    printf("End time = %f, N_step_total = %d, dt = %2.6e\n",
            End_Time, N_step_total, dt_host);
    
 #ifdef OPEN_MP
    itime = omp_get_wtime();
#else
    start_time = clock();
#endif
    // Compute the initial potential energy
    PE_total = Compute_PotentialEnergy_CPU(h_Pos, h_Mass, h_E_Potential, N_PARTICLE);
    // Compute the initial kinetic energy
    KE_total = Compute_KineticEnergy_CPU  (h_Vel, h_Mass, h_E_Kinetic  , N_PARTICLE);
    // Compute the total energy
    Energy_total = PE_total + KE_total;
    // Print the initial total energy
    printf("Initial Total Energy = %2.6e\n", Energy_total);
    // Save the initial snapshot
    SaveDataHDF5(h_Pos, h_Vel, h_Acc, h_Mass, h_E_Potential, h_E_Kinetic,
                 N_PARTICLE, num_round, dt_host, time_host);

    time_host = 0.0; // Make sure time_host starts from 0.0
    // Start the evolution loop
    for (uint round = 1; round <= N_step_total; round++)
    {
        num_round = round;
        // Compute acc(t_i) using Pos(t_i)
        Compute_ACC_CPU(h_Pos, h_Acc, h_Mass, N_PARTICLE);
        // Update Velocity to t_i+1/2 (Kick)
        UpdateVelocity_CPU(dt_host * const_half, h_Vel, h_Acc, N_PARTICLE);
        // Update Position to t_i+1 (Drift)
        UpdatePosition_CPU(dt_host, h_Pos, h_Vel, N_PARTICLE);
        // Compute acc(t_i+1) using Pos(t_i+1)
        Compute_ACC_CPU(h_Pos, h_Acc, h_Mass, N_PARTICLE);
        // Update Velocity to t_i+1 (Kick)
        UpdateVelocity_CPU(dt_host * const_half, h_Vel, h_Acc, N_PARTICLE);
        // Update simulation time
        time_host += dt_host;
        // check if the time has reached the end time
        if (time_host >= End_Time){
            printf("time_host = %2.6f, dt = %2.6e\n", time_host, dt_host);
            printf("Code time has reached the End_Time, stop now.\n");
            break;
        }
        // check if we save a snapshot in this round
        if (round % N_DUMP_STEP == 0){
            printf("Total round = %d, current round = %d\n", N_step_total, round);
            printf("Saving a snapshot at time_host = %2.6f, dt = %2.6e\n",
                    time_host, dt_host);
            // Compute the potential energy
            PE_total = Compute_PotentialEnergy_CPU(h_Pos, h_Mass, h_E_Potential, N_PARTICLE);
            // Compute the kinetic energy
            KE_total = Compute_KineticEnergy_CPU  (h_Vel, h_Mass, h_E_Kinetic  , N_PARTICLE);
            // Compute the total energy
            Energy_total = PE_total + KE_total;
            // Print the total energy
            printf("Total Energy = %2.6e\n", Energy_total);
            // Save the snapshot
            SaveDataHDF5(h_Pos, h_Vel, h_Acc, h_Mass, h_E_Potential, h_E_Kinetic,
                         N_PARTICLE, num_round, dt_host, time_host);
            printf("Snapshot saved successfully!\n");
        }
    }
    // Compute the potential energy
    PE_total = Compute_PotentialEnergy_CPU(h_Pos, h_Mass, h_E_Potential, N_PARTICLE);
    // Compute the kinetic energy
    KE_total = Compute_KineticEnergy_CPU  (h_Vel, h_Mass, h_E_Kinetic  , N_PARTICLE);
    // Compute the total energy
    Energy_total = PE_total + KE_total;
    // Print the total energy
    printf("Final Total Energy at t = %2.6e\n", Energy_total, time_host);
    SaveDataHDF5(h_Pos, h_Vel, h_Acc, h_Mass, h_E_Potential, h_E_Kinetic,
        N_PARTICLE, num_round, dt_host, time_host);
    printf("The final snapshot saved successfully!\n");

#ifdef OPEN_MP
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    printf("This simulation costs %f seconds with %d threads\n", exec_time, nProcessors);
#else
    end_time = clock();
    printf("This simulation costs %f seconds without OpenMP (single thread mode)\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

    // Release memory
    delete[] h_Mass;
    delete[] h_E_Potential;
    delete[] h_E_Kinetic;
    delete[] h_Pos;
    delete[] h_Vel;
    delete[] h_Acc;
    return EXIT_SUCCESS;        
}   // End of main function
