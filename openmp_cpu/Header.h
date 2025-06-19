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

#endif
