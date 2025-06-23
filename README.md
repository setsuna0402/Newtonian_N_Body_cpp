# Readme

## description:
This repository demonstrates a Newtonian N-body simulation with a softening length.
The codes are written in C++ and support shared-memory parallelisation by OpenMP and 
GPU concurrent parallelisation by CUDA.

Here are some assumptions adopted in these codes:
Use Leapfrog (kick-drift-kick) to update the position and velocity.
Use the softened gravity method to calculate the acceleration
The length of the timestep is a constant. (No adaptive timestep)
The Gravitational Constant is set to be 1, so the mass unit is determined by distance and time units.
Users need to make sure the unit conversion factors are correct:
G_Length_CodeUnits: convert distance from code unit to meter
G_Time_CodeUnits: convert time from code unit to second
G_Mass_CodeUnits: convert mass from code unit to kg
G_Velocity_CodeUnits: convert velocity from code unit to meter/second
G_Energy_CodeUnits: convert energy from code unit to Joule
For validating the dynamical algorithm, you can set all the code units to be one.
Coordinate order: [x, y, z]

For the OpenMP version (stored in the Folder openmp_cpu), you are free to change the following macros to configure your simulations:
#define FLOAT8  (on: run in double precision; otherwise run in single precision)
#define OPEN_MP (on: multiple threads (OpenMP))
#define SIMD    (on: use SIMD (AVX2, AVX512, etc. up to your cpu and compile configuration))
#define OPENMP_NUM_THREAD 8 (Set the Number of OpenMPthread)

For the Cuda version (stored in the Folder cuda_gpu), you are free to change the following macros to configure your simulations:
#define FLOAT8   (on: run in double precision; otherwise run in single precision, which is highly recommended for low-end GPUs)
#define OPEN_MP  (not in use in current version)
#define SIMD     (not in use in current version)
You need to choose ONE (and only one) of the following three options to choose different GPU implementations. 
#define GPU_SLOW   (test GPU code (simplest and slowest methods))
#define GPU_FAST   (use **per-thread** registers memory in GPU code (Faster))
#define GPU_SHARED (use **per-thread** registers and shared memory in GPU code (Fastest))
#define GPU_BLOCK_SIZE  (Number of threads in ONE block)

Current status:
Both versions demonstrate energy conservation (numerically) in test problems.

Future:
Currently, the GPU version supports only one main stream. It is worth trying to make it support multiple streams with OpenMP.





