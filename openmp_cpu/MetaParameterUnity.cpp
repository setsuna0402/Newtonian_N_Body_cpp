/***********************************************************************
/
/  Set all the code units to be one
/  Since the gravitational constant is assumed to be 1,
/  setting all code units to one makes the magnitude of gravity not physically meaningful,
/  but it simplifies the calculations and comparisons.
/
/  written by: KH 
/  date:       2025/6/19

************************************************************************/

#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Set all the code units to be one
int MetaParameterUnity(){

    cout << "Which function is running?     " << __func__ << endl;
#ifdef FLOAT8
    G_Length_CodeUnits   = 1.0;
    G_Time_CodeUnits     = 1.0;
    G_Mass_CodeUnits     = 1.0;
    G_Velocity_CodeUnits = 1.0;
    G_Energy_CodeUnits   = 1.0;
#else
    G_Length_CodeUnits   = 1.0f;    
    G_Time_CodeUnits     = 1.0f;
    G_Mass_CodeUnits     = 1.0f;
    G_Velocity_CodeUnits = 1.0f;
    G_Energy_CodeUnits   = 1.0f;
#endif
    // Print the conversion factors
    cout << "All code units are set to one!" << endl;
    cout << "G_Length_CodeUnits   = " << G_Length_CodeUnits   << endl;
    cout << "G_Time_CodeUnits     = " << G_Time_CodeUnits     << endl;
    cout << "G_Mass_CodeUnits     = " << G_Mass_CodeUnits     << endl;
    cout << "G_Velocity_CodeUnits = " << G_Velocity_CodeUnits << endl;
    cout << "G_Energy_CodeUnits   = " << G_Energy_CodeUnits   << endl;
    cout << __func__ << "...done!" << endl;
    // Return 0 to indicate success
    return 0;
}
