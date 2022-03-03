#include "utils.h"

#ifndef USE_MPI
/* compute kinetic energy */
void ekin(mdsys_t *sys,particles_t *part)
{
    double ekin = 0.0;
    // parallelization of outer loop
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:ekin)
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        ekin += 0.5*mvsq2e*sys->mass*(part[i].vx*part[i].vx + part[i].vy*part[i].vy + part[i].vz*part[i].vz);
    }

    sys->ekin = ekin;
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
#else
void ekin(mdsys_t *sys, particles_t *part)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<myq; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(part[i].vx*part[i].vx + part[i].vy*part[i].vy + part[i].vz*part[i].vz);
    }

    MPI_Allreduce(MPI_IN_PLACE,&sys->ekin,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // printf("Kinetic Energy %lf\n",sys->ekin);
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;

}
#endif