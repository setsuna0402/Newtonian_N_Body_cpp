/***********************************************************************
/
/  Compute the timestep which is safe for the simulation
/  In a two body system (with identical particle mass), 
/  the maximum acceleration is given by:
/  a_max = 2 / (epsion^2 * 3 ^(3/2))
/  delta_v = a_max * delta_t
/  The time step is defined by:
/  delta_t_safe = TimeSafetFactor * 2.6 * epsion ^ 2.0
/  where TimeSafetFactor is a safety factor to ensure the stability of the simulation
/  TimeSafetFactor ~ 0.01 is OK for a two body system (with identical particle mass)
/  Depending on the system, you may adopt a different way to compute the time step
/  written by: KH 
/  date:       2025/6/20
************************************************************************/
#include "macros_and_parameters.h"
#include "Header.h"

using namespace std;

// Compute TimeStep based on the theoretical maximum acceleration
void Safest_TimeStep_CPU(real &dt){
    static uint call_counter = 0;
    if (call_counter == 0) cout << "Which function is running?   " << __func__ << endl;
    const real eps2 = SOFTEN * SOFTEN; // the soften term in the soften gravity method
    dt = TimeSafety * 2.6 * eps2;
    printf("dt = %e (code unit)\n", dt);
    if (call_counter == 0) cout << __func__ << "...done!" << endl;
    call_counter++;
} // FUNCTION : Safest_TimeStep_CPU
