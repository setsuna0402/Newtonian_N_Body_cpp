/***********************************************************************
/
/  Load initial condition (HDF5 format)
/  Suggestion: Use double precision for the initial condition
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Mass (Dynamical array) : Mass of each particle (Size: N)
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Read the initial condition (HDF5 format)
int ICDataReaderHDF5(real (*Pos)[3], real (*Vel)[3], real *Mass){
    cout << "Which function is running?     " << __func__ << endl;
    int    rank          = 0;
    int    status_n      = 0;
    hid_t   file, dataset, group, attr;
    hid_t   datatype, dataspace       ;
    hsize_t dimsf_vector [2];  // dataset dimensions (2D array)
    hsize_t dimsf_scaler [1];  // dataset dimensions (1D array)
    herr_t  status ,ret     ;  // Return value
    //    size_t  size                   //dataspace handle
    //open the hdf5 file
    file  = H5Fopen(IC_FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        printf("Failed to open file {%s}\n", IC_FILE_NAME);
        printf("Abort\n");
        return EXIT_FAILURE;
    }
    // read Position
    dataset   = H5Dopen2(file, DATASETNAME_Pos, H5P_DEFAULT)    ;
    dataspace = H5Dget_space(dataset);    //dataspace handle
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dimsf_vector, NULL);
    cout<<"Position:  rank = "<<rank<<"  dimension  = "<<dimsf_vector[0]<<"x"<<dimsf_vector[1]<<endl;
    if (dimsf_vector[0] != N_PARTICLE || dimsf_vector[1] != 3) {
        printf("The dimension of Position is not correct!\n");
        printf("It should be %lld x 3, but it is %d x %d\n", N_PARTICLE, dimsf_vector[0], dimsf_vector[1]);
        printf("Abort\n");
        return EXIT_FAILURE;
    }
#ifdef FLOAT8
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, Pos);
#else
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, Pos);
#endif
    if (status != 0) {
        printf("H5Dread status = %d, it doesn't equal to zeros \n", status);
        printf("Failed to read {%s}\n", DATASETNAME_Pos);
        printf("Abort\n");
        return EXIT_FAILURE;
    }
    // read Velocity
    dataset   = H5Dopen2(file, DATASETNAME_Pos, H5P_DEFAULT)    ;
    dataspace = H5Dget_space(dataset);    //dataspace handle
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dimsf_vector, NULL);
    cout<<"TotalEnergy:  rank = "<<rank<<"  dimension  = "<<dimsf_vector[0]<<"x"<<dimsf_vector[1]<<endl;
    if (dimsf_vector[0] != N_PARTICLE || dimsf_vector[1] != 3) {
        printf("The dimension of Velocity is not correct!\n");
        printf("It should be %lld x 3, but it is %d x %d\n", N_PARTICLE, dimsf_vector[0], dimsf_vector[1]);
        printf("Abort\n");
        return EXIT_FAILURE;
    }
#ifdef FLOAT8
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, Vel);
#else
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, Vel);
#endif
    if (status != 0) {
        printf("H5Dread status = %d, it doesn't equal to zeros \n", status);
        printf("Failed to read {%s}\n", DATASETNAME_Pos);
        printf("Abort\n");
        return EXIT_FAILURE;
    }
    // read Mass
    dataset   = H5Dopen2(file, DATASETNAME_Mass, H5P_DEFAULT)   ;
    dataspace = H5Dget_space(dataset);    //dataspace handle
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dimsf_scaler, NULL);
    cout<<"Mass:  rank = "<<rank<<"  dimension  = "<<dimsf_scaler[0]<<endl;
    if (dimsf_scaler[0] != N_PARTICLE) {
        printf("The dimension of Mass is not correct!\n");
        printf("It should be %lld, but it is %d\n", N_PARTICLE, dimsf_scaler[0]);
        printf("Abort\n");
        return EXIT_FAILURE;
    }   
#ifdef FLOAT8
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, Mass);
#else
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, Mass);
#endif
    if (status != 0) {
        printf("H5Dread status = %d, it doesn't equal to zeros \n", status);
        printf("Failed to read {%s}\n", DATASETNAME_Mass);
        printf("Abort\n");
        return EXIT_FAILURE;
    }
                     
    // Close the dataset and file
    // H5Aclose(attr   );
    // H5Gclose(group  );
    H5Dclose(dataset);
    H5Fclose(file   );
    cout << __func__ << "...done!" << endl;
    return 0;
} // ICDataReaderHDF5