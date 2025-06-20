/***********************************************************************
/
/  Create a random initial condition for particles in a 3D space
/
/  written by: KH 
/  date:       2025/6/20
/  Input:
/  Max_R (const real) : Maximum radius of the particles (code unit)
/  Pos  (Dynamical array) : Position of each particle (Size: N x 3)
/  Vel  (Dynamical array) : Velocity of each particle (Size: N x 3)
/  Size (const uint) : Number of particles
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Create a random initial condition for particles in a 3D space
void CreateRandomIC(const real Max_R, real (*Pos)[3], real (*Vel)[3], const uint Size){
    static uint call_counter = 0;
    if (call_counter == 0) cout << "Which function is running?   " << __func__ << endl;
    real Radius = 0.0;
    const real centre[3] = {0.5 * BOXSIZE_CODE, 0.5 * BOXSIZE_CODE, 0.5 * BOXSIZE_CODE};
    /*the time cost here is so low, parallel this loop will solw down the code..*/
    /*
    #ifdef OPEN_MP
    #pragma omp parallel private(Radius)
    {
        #pragma omp for private(Radius)
    #endif
    */
    for (int i = 0; i < Size; i++){
        Radius = Max_R + 1.0;
        //    ensure r <= MaxR
        while (Radius > Max_R){
            Radius = 0.0;
            for (int d = 0; d < 3; d++){
                Pos[i][d] = ((real)rand() / (real)RAND_MAX) * 
                            2.0 * Max_R - Max_R + centre[d];
            }
        #ifdef FLOAT8
            for (int d = 0; d < 3; d++) Radius += pow((Pos[i][d] - centre[d]), 2.0);
            Radius = sqrt(Radius);
        #else
            for (int d = 0; d < 3; d++) Radius += powf((Pos[i][d] - centre[d]), 2.0f);
            Radius = sqrtf(Radius);
        #endif
        }
        //    -V_Max < v < V_Max
        for (int d = 0; d < 3; d++){
            Vel[i][d] = ((real)rand() / (real)RAND_MAX) * 2.0 * V_MAX - V_MAX;
        }
    }
    /*
    #ifdef OPEN_MP
    }
    #endif
    */
    if (call_counter == 0) cout << __func__ << "...done!" << endl;
    call_counter++;
}// FUNCTION : CreateRandomIC