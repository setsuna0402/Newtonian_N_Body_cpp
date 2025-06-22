/***********************************************************************
/  Set GPU globle variables  
/
/  written by: KH 
/  date:       2025/6/20
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Gobal GPU constants
__constant__ real d_Gc_GravConst_mks;  // Gravitational constant in mks unit
__constant__ real d_Gc_GravConst_cgs;  // Gravitational constant in cgs unit
__constant__ real d_Gc_sqrt_pi      ;  // sqrt(pi)
__constant__ real d_Gc_Mpc_to_km    ;  // 1Mpc = 3.08567758128E+19 km
__constant__ real d_Gc_Mpc_to_m     ;  // 1Mpc = 3.08567758128E+22 m
__constant__ real d_Gc_Mpc_to_cm    ;  // 1Mpc = 3.08567758128E+24 cm
__constant__ real d_Gc_Boxsize      ;  // in code unit
__constant__ real d_Gc_End_time     ;  // Endtime of simulation (code unit)
// Gobal GPU variables
__device__ real d_G_Length_CodeUnits    ;  // Convert distance from code unit to meter
__device__ real d_G_Time_CodeUnits      ;  // Convert time from code unit to second
__device__ real d_G_Mass_CodeUnits      ;  // Convert mass from code unit to kg
__device__ real d_G_Velocity_CodeUnits  ;  // Convert velocity from code unit to meter/second
__device__ real d_G_Energy_CodeUnits    ;  // Convert energy from code unit to Joule

// Pass globle variables from CPU to GPU
int PassGlobleVariabletoDevice(){
    cout << "Which function is running?     " << __func__ << endl;
    // __constant__ globle variables in GPU
    cudaMemcpyToSymbol(d_Gc_GravConst_mks , &Gc_GravConst_mks , sizeof(real));
    cudaMemcpyToSymbol(d_Gc_GravConst_cgs , &Gc_GravConst_cgs , sizeof(real));
    cudaMemcpyToSymbol(d_Gc_sqrt_pi       , &Gc_sqrt_pi       , sizeof(real));
    cudaMemcpyToSymbol(d_Gc_sqrt_pi       , &Gc_sqrt_pi       , sizeof(real));  
    cudaMemcpyToSymbol(d_Gc_Mpc_to_km     , &Gc_Mpc_to_km     , sizeof(real));  
    cudaMemcpyToSymbol(d_Gc_Mpc_to_m      , &Gc_Mpc_to_m      , sizeof(real));  
    cudaMemcpyToSymbol(d_Gc_Mpc_to_cm     , &Gc_Mpc_to_cm     , sizeof(real));  
    cudaMemcpyToSymbol(d_Gc_Boxsize       , &Gc_Boxsize       , sizeof(real)); 
    cudaMemcpyToSymbol(d_Gc_End_time      , &Gc_End_time      , sizeof(real)); 
    // __device__ globle variable in GPU 
    cudaMemcpyToSymbol(d_G_Length_CodeUnits   , &G_Length_CodeUnits   , sizeof(real)); 
    cudaMemcpyToSymbol(d_G_Time_CodeUnits     , &G_Time_CodeUnits     , sizeof(real)); 
    cudaMemcpyToSymbol(d_G_Mass_CodeUnits     , &G_Mass_CodeUnits     , sizeof(real)); 
    cudaMemcpyToSymbol(d_G_Velocity_CodeUnits , &G_Velocity_CodeUnits , sizeof(real));
    cudaMemcpyToSymbol(d_G_Energy_CodeUnits   , &G_Energy_CodeUnits   , sizeof(real));
    cudaDeviceSynchronize();
    cout << __func__ << "...done!" << endl;
    return 0;
} // FUNCTION : PassGlobleVariabletoDevice


