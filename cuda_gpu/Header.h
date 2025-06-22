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
// #include "/home/edward/Library/hdf5/1.10.3/include/hdf5.h"
#include "hdf5.h"
#include "macros_and_parameters.h"

// Read the meta parameter file. Mainly used to set the code units
int MetaParameterReader();
// Set all the code units to be one
int MetaParameterUnity();
// Compute the acceleration for each particles by soften gravity
void Compute_ACC_CPU(real (*Pos)[3], real (*Acc)[3], real *Mass, const uint Size);
// Compute TimeStep based on the theoretical maximum acceleration
void Safest_TimeStep_CPU(real &dt);
// Update position by a given timestep
void UpdatePosition_CPU(real dt, real (*Pos)[3], real (*Vel)[3], const uint Size);
// Update velocity by a given timestep
void UpdateVelocity_CPU(real dt, real (*Vel)[3], real (*Acc)[3], const uint Size);
// Create a random initial condition for particles in a 3D space
void CreateRandomIC(const real Max_R, real (*Pos)[3], real (*Vel)[3],
                    real *Mass, const uint Size);
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
// Save snapshot (HDF5 format)
void SaveDataHDF5(real (*Pos)[3], real (*Vel)[3], real (*Acc)[3], real *Mass,
                  real *E_Potential, real *E_Kinetic, const uint Size,
                  const uint Num_step, const real TimeStep, const real Time);
// Read the initial condition (HDF5 format)
int ICDataReaderHDF5(real (*Pos)[3], real (*Vel)[3], real (*Acc)[3], real *Mass);

#endif
