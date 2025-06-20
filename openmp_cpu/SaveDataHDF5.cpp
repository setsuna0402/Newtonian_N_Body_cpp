/***********************************************************************
/
/  Save data in HDF5 format
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Mass (Dynamical array) : Mass of each particle (Size: N)
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
/  Acc  (Dynamical array) : Acceleration of each particle (Size: N x 3)
/  E_Potential (Dynamical array) : Potential energy of each particle (Size: N)
/  E_Kinetic (Dynamical array) : Kinetic energy of each particle (Size: N)
/  Size (const uint ) : Number of particles
/  Num_step (const uint) : The number of evolved steps
/  TimeStep (const real) : The time step of the simulation
/  Time (const real) : The total time of the simulation
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

void SaveDataHDF5(real (*Pos)[3], real (*Vel)[3], real (*Acc)[3], real *Mass,
                  real *E_Potential, real *E_Kinetic, const uint Size,
                  const uint Num_step, const real TimeStep, const real Time){
    cout << "Which function is running?     " << __func__ << endl;
    const int RANK = 2;
    hid_t file, dataset;       // file and dataset handles
    hid_t datatype, dataspace; // handles
    hsize_t dimsf[RANK];       // dataset dimensions (2D array)
    hsize_t dimsf_time[1];     // dataset dimensions (1D array)
    hsize_t dimsf_scaler[1];   // dataset dimensions (1D array)
    herr_t status;
    string str1(OUTPUT_FILE_Prefix);
    string str2 = to_string(Num_step);
    string str3(".hdf5");
    str1.append(str2);
    str1.append(str3);
    cout << str1 << endl;
    cout << str2 << endl;
    file = H5Fcreate(str1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // Dimensions for the 2D datasets (Pos, Vel, Acc)
    dimsf[0] = Size;
    dimsf[1] = 3;
    // Dimensions for the 1D datasets (Mass, E_Potential, E_Kinetic)
    dimsf_scaler[0] = Size;
    // Dimensions for the 1D dataset (Time)
    dimsf_time[0] = 1;
#ifdef FLOAT8
    dataspace = H5Screate_simple(RANK, dimsf, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    // Save Positon
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Pos, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Pos);
    // Save Velocity
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Vel, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Vel);
    // Save Acceleration
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Acc, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Acc);
    // Save Mass
    dataspace = H5Screate_simple(RANK, dimsf_scaler, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Mass, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Mass);
    // Save Potential Energy
    dataspace = H5Screate_simple(RANK, dimsf_scaler, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_PE, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, E_Potential);
    // Save Kinetic Energy
    dataspace = H5Screate_simple(RANK, dimsf_scaler, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_KE, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, E_Kinetic);
    // Save TimeStep
    dataspace = H5Screate_simple(1, dimsf_time, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, "delta_t", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &TimeStep);
    // Save Simulation Time (code time)
    dataspace = H5Screate_simple(1, dimsf_time, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, "code_time", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Time);
#else
    dataspace = H5Screate_simple(RANK, dimsf, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_FLOAT);
    // Save Position
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Pos, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Pos);
    // Save Velocity
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Vel, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Vel);
    // Save Acceleration
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Acc, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Acc);
    // Save Mass
    dataspace = H5Screate_simple(RANK, dimsf_scaler, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_FLOAT);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_Mass, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Mass);
    // Save Potential Energy
    dataspace = H5Screate_simple(RANK, dimsf_scaler, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_FLOAT);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_PE, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, E_Potential);
    // Save Kinetic Energy
    dataspace = H5Screate_simple(RANK, dimsf_scaler, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_FLOAT);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, DATASETNAME_KE, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, E_Kinetic);
    // Save TimeStep
    dataspace = H5Screate_simple(1, dimsf_time, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_FLOAT);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, "delta_t", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &TimeStep);
    // Save Simulation Time (code time)
    dataspace = H5Screate_simple(1, dimsf_time, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_FLOAT);
    status    = H5Tset_order(datatype, H5T_ORDER_LE);
    dataset   = H5Dcreate2(file, "code_time", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Time);
#endif

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    cout << __func__ << "...done!" << endl;
} // Function      Writer_potential_HDF5