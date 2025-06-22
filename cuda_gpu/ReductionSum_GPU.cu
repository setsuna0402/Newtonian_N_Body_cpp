/***********************************************************************
/
/  Reduction sum over the target array  (GPU version)
/  The return values (an array) are the local sum of each GPU Block
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Target (Dynamical array) : 1D array  (Size: N)
/  Output (Dynamical array) : The local summation of each GPU block (Size: (int)N/Block_size + 1)
/  Size (const uint ) : The number of elements in the Target array

************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"


// GPU function : Reduction sum over the target array
__global__ void ReductionSum_GPU(real *Target, real *Output, const uint Size){
    // the index of the thread in GPU
    const int idx = blockDim.x * blockIdx.x + threadIdx.x;
    // the index of the thread in the block
    const int tid = threadIdx.x;
    // shared memory for each block
    __shared__ real sdata[GPU_BLOCK_SIZE];
    
    // load data into shared memory
    __syncthreads(); // synchronise all threads before loading data
    if (idx < Size) {
        sdata[tid] = Target[idx];
    } else {
        sdata[tid] = 0.0f; // avoid if-else branching in the loop
    }
    __syncthreads(); // ensure all threads have loaded their data

    // perform reduction in shared memory
    // Assuming that GPU_BLOCK_SIZE is a power of 2
    for (int stride = GPU_BLOCK_SIZE / 2; stride > 0; stride /= 2) {
        if (tid < stride) {
            sdata[tid] += sdata[tid + stride];
        }
        __syncthreads(); // synchronize threads after each reduction step
    }
    // write the result of this block to the output array
    if (tid == 0) {
        Output[blockIdx.x] = sdata[0];
    }
} // Function: ReductionSum_GPU
