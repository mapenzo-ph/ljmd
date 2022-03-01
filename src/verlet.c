#include "utils.h"
/* velocity verlet */
void verlet1(mdsys_t *sys)
{
    /* first part: propagate velocities by half and positions by full step */

    // add openMP parallel loop
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }
}

void verlet2(mdsys_t *sys)
{
    /* second part: propagate velocities by another half step */

    // add openMP parallel loop
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}
