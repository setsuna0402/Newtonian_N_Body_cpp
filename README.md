This repository demostrates Newtonian N body simulation with soften length.
The codes write in C++ and supports shared-memory parallelisation by OpenMP and 
gpu concurrent parallelisation by CUDA.

Here are some assumptions adopted in these codes:
Use Leapfrog (kick-drift-kick) to update the position and velocity.
Use soften gravity method to calculate the acceleration
Length of timestep is a constant. (No adaptive timestep)
Gravitational Constant is set to be 1, so mass unit is determinated by distance and time units.
User need to make sure the unit conversion factors are correct:
G_Length_CodeUnits: convert distance from code unit to meter
G_Time_CodeUnits  : convert time from code unit to second
G_Mass_CodeUnits  : convert mass from code unit to kg
G_Velocity_CodeUnits: convert velocity from code unit to meter/second
G_Energy_CodeUnits: convert energy from code unit to Joule
For validating the dynamical algorithm, you can set all the code units to be one.
Coordinate order: [x, y, z]

For openmp version (stored in Folder openmp_cpu), you are free to change the following macros to configure youre simulations:
#define FLOAT8  (on: run in double precision; otherwise run in single precision)
#define OPEN_MP (on: multiple threads (OpenMP))
#define SIMD    (on: use SIMD (AVX2, AVX512, etc. up to your cpu and compile configuration))
#define OPENMP_NUM_THREAD 8 (Set the Number of OpenMPthread)

For cuda version (stored in Folder cuda_gpu), you are free to change the following macros to configure youre simulations:
#define FLOAT8   (on: run in double precision; otherwise run in single precision, which is highly recommended for low-end GPUs)
#define OPEN_MP  (not in used in current version)
#define SIMD     (not in used in current version)
You need to choose ONE (and only one) of the following three options to choose different GPU implementations. 
#define GPU_SLOW   (test GPU code (simplest and slowest methods))
#define GPU_FAST   (use **per-thread** registers memory in GPU code (Faster))
#define GPU_SHARED (use **per-thread** registers and shared memory in GPU code (Fastest))
#define GPU_BLOCK_SIZE  (Number of theard in ONE block)

Current status:
Both versions demostrate energy conservation (numerically) in test problems.

Future:
Currenty, GPU version supports only one main stream. It is worth to try to make it support multiple streams with OpenMP.





