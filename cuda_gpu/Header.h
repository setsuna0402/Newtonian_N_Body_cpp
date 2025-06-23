#ifndef HEADER_H
#define HEADER_H
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <omp.h>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
// #include "/home/edward/Library/hdf5/1.10.3/include/hdf5.h"
#include "hdf5.h"
#include "macros_and_parameters.h"

// Read the meta parameter file. Mainly used to set the code units
int MetaParameterReader();
// Set all the code units to be one
int MetaParameterUnity();
// Compute TimeStep based on the theoretical maximum acceleration
void Safest_TimeStep_CPU(real &dt);
// Create a random initial condition for particles in a 3D space
void CreateRandomIC(const real Max_R, real (*Pos)[3], real (*Vel)[3],
                    real *Mass, const uint Size);
// Save snapshot (HDF5 format)
void SaveDataHDF5(real (*Pos)[3], real (*Vel)[3], real (*Acc)[3], real *Mass,
                  real *E_Potential, real *E_Kinetic, const uint Size,
                  const uint Num_step, const real TimeStep, const real Time);
// Read the initial condition (HDF5 format)
int ICDataReaderHDF5(real (*Pos)[3], real (*Vel)[3], real (*Acc)[3], real *Mass);
// Pass globle variables from CPU to GPU
int PassGlobleVariabletoDevice();
// Function to check CUDA errors
void checkCuda(cudaError_t result, const char *msg = "",
               const char *file = __FILE__, int line = __LINE__);
// GPU implementation functions
// GPU function : Compute the acceleration for each particles
__global__ void Compute_ACC_GPU_SLOW(real (*Pos)[3], real (*Acc)[3],
                                     real *Mass, const uint Size);
// GPU function : Compute the acceleration (using per-thread register) for each particles
__global__ void Compute_ACC_GPU_FAST(real (*Pos)[3], real (*Acc)[3],
                                     real *Mass, const uint Size);
// GPU function : Compute the acceleration (using shared memory) for each particles
__global__ void Compute_ACC_GPU_SHARED(real (*Pos)[3], real (*Acc)[3],
                                       real *Mass, const uint Size);
// GPU function : Update the position of all the particles by a given timestep
__global__ void UpdatePosition_GPU(real dt, real (*Pos)[3], real (*Vel)[3],
                                   const uint Size);
// GPU function : Update the Velocity of all the particles by a given timestep
__global__ void UpdateVelocity_GPU(real dt, real (*Vel)[3], real (*Acc)[3],
                                   const uint Size);
// GPU function : Compute the potential for each particles
__global__ void Compute_PotentialEnergy_GPU_SLOW(real (*Pos)[3], real *Mass,
                                                 real *E_Potential, const uint Size);
// GPU function : Compute the potential for each particles (using per-thread register)
__global__ void Compute_PotentialEnergy_GPU_FAST(real (*Pos)[3], real *Mass,
                                                 real *E_Potential, const uint Size);
// GPU function : Compute the potential for each particles (using shared memory)
__global__ void Compute_PotentialEnergy_GPU_SHARED(real (*Pos)[3], real *Mass,
                                                   real *E_Potential, const uint Size);
// GPU function : Compute the kinetic energy for each particles
__global__ void Compute_KineticEnergy_GPU(real (*Vel)[3], real *Mass,
                                          real *E_Kinetic, const uint Size);  
// GPU function : Reduction sum over the target array
__global__ void ReductionSum_GPU(real *Target, real *Output, const uint Size);                                          
                                                                                           

/*
// Some functions for CPU implementation
// Compute the acceleration for each particles by soften gravity
void Compute_ACC_CPU(real (*Pos)[3], real (*Acc)[3], real *Mass, const uint Size);
// Update position by a given timestep
void UpdatePosition_CPU(real dt, real (*Pos)[3], real (*Vel)[3], const uint Size);
// Update velocity by a given timestep
void UpdateVelocity_CPU(real dt, real (*Vel)[3], real (*Acc)[3], const uint Size);
// Compute the potential energy for each particles
real Compute_PotentialEnergy_CPU(real (*Pos)[3], real *Mass,
                                 real *E_Potential, const uint Size);
// Compute the total potential energy of the system
real ComputeTotalPotentialEnergy_CPU(real *E_Potential, const uint Size);  
// Compute the kinetic energy for each particles by soften gravity
real Compute_KineticEnergy_CPU(real (*Vel)[3], real *Mass, real *E_Kinetic,
                               const uint Size);
// Compute the total kinetic energy of the system
real ComputeTotalKineticEnergy_CPU(real *E_Kinetic, const uint Size);
*/

#endif
