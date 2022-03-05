#include "utils.h"

void verlet1(mdsys_t *sys, coords_t * coord, vel_t *vel, for_t* forces)
{
    int i;
    double tmp=0.5*sys->dt/mvsq2e;

    /* first part: propagate velocities by half and positions by full step */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (i=0; i<sys->natoms; ++i) {
        // Update velocities
        vel[i].vx += tmp * forces[i].fx/sys->mass;
        vel[i].vy += tmp * forces[i].fy/sys->mass;
        vel[i].vz += tmp * forces[i].fz/sys->mass;
    }
    for(i=0; i<sys->natoms;++i){
        // Update positions
        coord[i].rx += sys->dt*vel[i].vx;
        coord[i].ry += sys->dt*vel[i].vy;
        coord[i].rz += sys->dt*vel[i].vz;
    }
}


void verlet2(mdsys_t *sys, vel_t *vel, for_t* forces){
    double tmp=0.5*sys->dt/mvsq2e;
    /* compute forces and potential energy */
    /* second part: propagate velocities by another half step */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<sys->natoms; ++i) {
        vel[i].vx += tmp * forces[i].fx/sys->mass;
        vel[i].vy += tmp * forces[i].fy/sys->mass;
        vel[i].vz += tmp * forces[i].fz/sys->mass;
    }
}
