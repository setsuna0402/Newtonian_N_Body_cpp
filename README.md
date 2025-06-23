# Readme

## Description:
This repository demonstrates a **Newtonian N-body** simulation with a **softening length**.
The codes are written in C++ and support shared-memory parallelisation by OpenMP and 
GPU concurrent parallelisation by CUDA.

Here are some assumptions adopted in these codes:

Use Leapfrog (kick-drift-kick) to update the position and velocity.

Use the softened gravity method to calculate the acceleration

The length of the timestep is a constant. (No adaptive timestep)

The Gravitational Constant is set to be 1, so the mass unit is determined by distance and time units.
Users have to make the unit conversion factors consistent. For validating the dynamical algorithm, you can set all the code units to be one.

## For the OpenMP version (stored in the Folder openmp_cpu): 
You are free to change the following macros to configure your simulations:

FLOAT8  (on: run in double precision; otherwise run in single precision)

OPEN_MP (on: multiple threads (OpenMP))

SIMD    (on: use SIMD (AVX2, AVX512, etc. up to your cpu and compile configuration))

OPENMP_NUM_THREAD 8 (Set the Number of OpenMPthread)

## For the Cuda version (stored in the Folder cuda_gpu):
You are free to change the following macros to configure your simulations:

FLOAT8   (on: run in double precision; otherwise run in single precision, which is highly recommended for low-end GPUs)

OPEN_MP  (not in use in the current version)

SIMD     (not in use in the current version)

You need to choose **ONLY ONE** of the following three options to adopt different GPU implementations. 

GPU_SLOW   (use the simplest and slowest methods for testing)

GPU_FAST   (use **per-thread** registers memory in GPU code (Faster))

GPU_SHARED (use **per-thread** registers and shared memory in GPU code (Fastest))

GPU_BLOCK_SIZE  256 (Number of threads in ONE block)

## Current status:
Both versions demonstrate energy conservation (numerically) in test problems.

## Future:
Currently, the GPU version supports only one main stream. It is worth trying to make it support multiple streams with OpenMP.





