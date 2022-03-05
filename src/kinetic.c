#include "utils.h"

/* compute kinetic energy */
void ekin(mdsys_t *sys,vel_t *vel)
{
    double ekin = 0.0;
    double tmp= 0.5*mvsq2e*sys->mass;
    // parallelization of outer loop
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:ekin)
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        ekin += tmp*(vel[i].vx*vel[i].vx + vel[i].vy*vel[i].vy + vel[i].vz*vel[i].vz);
    }

    sys->ekin = ekin;
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
