/***********************************************************************
/
/  Read Meta parameter file and set all the relative parameters
/  Set the code units for the simulation
/
/  written by: KH 
/  date:       2025/6/19

************************************************************************/

#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Read the meta parameter file. Mainly used to set the code units
int MetaParameterReader(){

    cout << "Which function is running?     " << __func__ << endl;
    unsigned int N_SIZE = 0;
    double temp = 0.0;
    hid_t file, group, attr;
    herr_t ret; // Return value
    file = H5Fopen(META_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    group = H5Gopen1(file, "Header");
    attr  = H5Aopen(group, "Number_of_Particles", H5P_DEFAULT ) ;
    ret   = H5Aread(attr , H5T_STD_I32LE   , &N_SIZE                 ) ;
    if (N_SIZE != N_PARTICLE) {
        fprintf(stderr, "The numebr of particle in metafile, N_SIZE = %d \n", N_SIZE);
        fprintf(stderr, "The numebr of particle in code, N_PARTICLE = %d \n", N_PARTICLE);
        fprintf(stderr, "N_PARTICLE doesn't equal to N_SIZE\n");
        return EXIT_FAILURE;
    }

#ifdef FLOAT8
    attr  = H5Aopen(group, "Length_CodeUnits"           , H5P_DEFAULT ) ;
    ret   = H5Aread(attr , H5T_IEEE_F64LE  , &G_Length_CodeUnits      ) ; 

    attr  = H5Aopen(group, "Time_CodeUnits"             , H5P_DEFAULT ) ;
    ret   = H5Aread(attr , H5T_IEEE_F64LE  , &G_Time_CodeUnits        ) ;
    // Calculate the mass code unit, assuming gravitational constant is 1
    G_Mass_CodeUnits     = pow(G_Length_CodeUnits, 3.0) / 
                           (pow(G_Time_CodeUnits, 2.0) * Gc_GravConst_mks);
    G_Velocity_CodeUnits = G_Length_CodeUnits / G_Time_CodeUnits ;
    G_Energy_CodeUnits   = G_Mass_CodeUnits * pow(G_Velocity_CodeUnits, 2.0);
#else
    // For single precision
    attr  = H5Aopen(group, "Length_CodeUnits"           , H5P_DEFAULT ) ;
    ret   = H5Aread(attr , H5T_IEEE_F32LE  , &G_Length_CodeUnits      ) ; 

    attr  = H5Aopen(group, "Time_CodeUnits"             , H5P_DEFAULT ) ;
    ret   = H5Aread(attr , H5T_IEEE_F32LE  , &G_Time_CodeUnits        ) ;

    // Calculate the mass code unit in double precision
    temp  = pow(G_Length_CodeUnits, 3.0) / (pow(G_Time_CodeUnits, 2.0) * Gc_GravConst_mks);
    G_Mass_CodeUnits     = (float)(temp); // Convert to single precision
    G_Velocity_CodeUnits = G_Length_CodeUnits / G_Time_CodeUnits ;
    G_Energy_CodeUnits   = G_Mass_CodeUnits * powf(G_Velocity_CodeUnits, 2.0);
                      
#endif
    H5Aclose(attr);
    H5Gclose(group);
    H5Fclose(file);
    // Print the conversion factors
    cout << "G_Length_CodeUnits   = " << G_Length_CodeUnits   << endl;
    cout << "G_Time_CodeUnits     = " << G_Time_CodeUnits     << endl;
    cout << "G_Mass_CodeUnits     = " << G_Mass_CodeUnits     << endl;
    cout << "G_Velocity_CodeUnits = " << G_Velocity_CodeUnits << endl;
    cout << "G_Energy_CodeUnits   = " << G_Energy_CodeUnits   << endl;
    cout << __func__ << "...done!" << endl;
    return 0;
}

