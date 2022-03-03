#include "utils.h"

#ifndef USE_MPI
/* velocity verlet */
void verlet1(mdsys_t *sys,particles_t *part)
{
    /* first part: propagate velocities by half and positions by full step */

    // add openMP parallel loop
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        part[i].vx += 0.5*sys->dt / mvsq2e * part[i].fx / sys->mass;
        part[i].vy += 0.5*sys->dt / mvsq2e * part[i].fy / sys->mass;
        part[i].vz += 0.5*sys->dt / mvsq2e * part[i].fz / sys->mass;
        part[i].rx += sys->dt*part[i].vx;
        part[i].ry += sys->dt*part[i].vy;
        part[i].rz += sys->dt*part[i].vz;
    }
}

void verlet2(mdsys_t *sys,particles_t *part)
{
    /* second part: propagate velocities by another half step */

    // add openMP parallel loop
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        part[i].vx += 0.5*sys->dt / mvsq2e * part[i].fx / sys->mass;
        part[i].vy += 0.5*sys->dt / mvsq2e * part[i].fy / sys->mass;
        part[i].vz += 0.5*sys->dt / mvsq2e * part[i].fz / sys->mass;
    }
}
#else
extern int myq;

void verlet1(mdsys_t *sys, particles_t *part)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<myq; ++i) {
        // Update velocities
        part[i].vx += 0.5*sys->dt / mvsq2e * part[i].fx / sys->mass;
        part[i].vy += 0.5*sys->dt / mvsq2e * part[i].fy / sys->mass;
        part[i].vz += 0.5*sys->dt / mvsq2e * part[i].fz / sys->mass;
        // Update positions
        part[i].rx += sys->dt*part[i].vx;
        part[i].ry += sys->dt*part[i].vy;
        part[i].rz += sys->dt*part[i].vz;
    }
}

void verlet2(mdsys_t *sys, particles_t *part){
    int i;
    /* second part: propagate velocities by another half step */
    for (i=0; i<myq; ++i) {
        part[i].vx += 0.5*sys->dt / mvsq2e * part[i].fx / sys->mass;
        part[i].vy += 0.5*sys->dt / mvsq2e * part[i].fy / sys->mass;
        part[i].vz += 0.5*sys->dt / mvsq2e * part[i].fz / sys->mass;
    }
}

#endif