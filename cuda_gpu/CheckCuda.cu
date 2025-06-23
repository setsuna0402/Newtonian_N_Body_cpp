/***********************************************************************
/
/  Check whether GPU function runs successfully
/  written by: KH 
/  date:       2025/6/23
/  Input:
/  result (cudaError_t) : The status of the CUDA function
/  msg (const char *) : The message to be printed if an error occurs 
/                       (usually the function name)
/  file (const char *) : The file name where the error occurs
/  line (int) : The line number where the error occurs

************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

// Function to check CUDA errors
void checkCuda(cudaError_t result, const char *msg = "", const char *file = __FILE__,
               int line = __LINE__){
    if (result != cudaSuccess) {
        printf("CUDA Error: %s, at file %s on line %d\n", msg, file, line);
        printf("CUDA Error code: %d (%s)\n", result, cudaGetErrorString(result));
        fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(result));
        fprintf(stderr, "Error (%s) occurred in file %s at line %d\n", msg, file, line);
        exit(EXIT_FAILURE);
    }
}