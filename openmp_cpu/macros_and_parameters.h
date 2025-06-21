#ifndef __macros_and_parameters_h_
#define __macros_and_parameters_h_

// single/double precision
#define FLOAT8  // on: double precision (If you read IC, use double precision)
// Use openMP for parallel computing
#define OPEN_MP // on: multiple threads (OpenMP); off: single thread 
// #define SIMD    // on: use SIMD (AVX2, AVX512, etc.)
#define DEBUG   // on: debug mode and print debug information

#ifdef OPEN_MP
    #define OPENMP_NUM_THREAD 8
#else
    #define OPENMP_NUM_THREAD 1
#endif

#ifdef FLOAT8
typedef double real;
#else
typedef float real; 
#endif

#define IC_FILE_NAME "IC_N_1024.hdf5"  // input file name (IC)
#ifdef FLOAT8
    #define META_NAME "N_body_meta_double.h5"       // input meta data file name 
#else
    #define META_NAME "N_body_meta_single.h5"       // input meta data file name 
#endif
#define OUTPUT_FILE_Prefix "N_body_8_thread_SIMD_N_1024_step_" // output file name(HDF5)
// Dataset names
#define DATASETNAME_Pos  "Position"
#define DATASETNAME_Vel  "Velocity"
#define DATASETNAME_Acc  "Acceleration"
#define DATASETNAME_Mass "Mass"
#define DATASETNAME_PE   "PotentialEnergy"
#define DATASETNAME_KE   "KineticEnergy"

#ifdef FLOAT8
    #define BOXSIZE_CODE    1.0     // Boxsize in code unit 
    #define FLOAT_EPSILON   1.0e-10 // Small number as a threshold
    #define CourantNumber   0.5     // Safety Number satisfying Courant Condition
    // soften length for calculating the gravitational acceleration
    #define SOFTEN          0.01   // Soften Length (code unit)
    #define End_Time        0.1     // Endtime of simulation (code unit)
    #define TimeSafety      1.0e-2  // safety factor for computing timestep
    #define V_MAX           0.1    // Maximum initial velocity of particles (code unit)
#else
    #define BOXSIZE_CODE    1.0f    // Boxsize in code unit 
    #define FLOAT_EPSILON   1.0e-5f // Small number as a threshold
    #define CourantNumber   0.5f    // Safety Number satisfying Courant Condition
    // soften length for calculating the gravitational acceleration
    #define SOFTEN          0.01f  // Soften Length (code unit)
    #define End_Time        0.1f    // Endtime of simulation (code unit)
    #define TimeSafety      1.0e-2f // safety factor for computing timestep
    #define V_MAX           0.1f   // Maximum initial velocity of particles (code unit)
#endif

#define N_PARTICLE      1024             // Number of particles
// Maximum number of evolution steps
#define N_STEP_MAX      1000000          // Maximum number of evolution steps
#define N_DUMP_STEP     10000            // Number of dump steps
#define N_DUMP_MAX      (N_STEP_MAX / N_DUMP_STEP) // Maximum number of dump files

// Gobal constants
extern const real Gc_GravConst_mks;  // Gravitational constant in mks unit
extern const real Gc_GravConst_cgs;  // Gravitational constant in cgs unit
extern const real Gc_sqrt_pi      ;  // sqrt(pi)
extern const real Gc_Mpc_to_km    ;  // 1Mpc = 3.08567758128E+19 km
extern const real Gc_Mpc_to_m     ;  // 1Mpc = 3.08567758128E+22 m
extern const real Gc_Mpc_to_cm    ;  // 1Mpc = 3.08567758128E+24 cm
extern const real Gc_Boxsize      ;  // in code unit
extern const real Gc_End_time     ;  // Endtime of simulation (code unit)
extern real G_Length_CodeUnits    ;  // Convert distance from code unit to meter
extern real G_Time_CodeUnits      ;  // Convert time from code unit to second
extern real G_Mass_CodeUnits      ;  // Convert mass from code unit to kg
extern real G_Velocity_CodeUnits  ;  // Convert velocity from code unit to meter/second
extern real G_Energy_CodeUnits    ;  // Convert energy from code unit to Joule

#endif
