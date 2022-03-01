#include "utils.h"

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{
    double ekin = 0.0;
    // parallelization of outer loop
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:ekin)
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }

    sys->ekin = ekin;
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
