#include "utils.h"

// /* velocity verlet */
// void verlet1(mdsys_t *sys,particles_t *part)
// {
//     /* first part: propagate velocities by half and positions by full step */

//     // add openMP parallel loop
//     #ifdef _OPENMP
//     #pragma omp parallel for
//     #endif
//     for (int i=0; i<sys->natoms; ++i) {
//         part[i].vx += 0.5*sys->dt / mvsq2e * part[i].fx / sys->mass;
//         part[i].vy += 0.5*sys->dt / mvsq2e * part[i].fy / sys->mass;
//         part[i].vz += 0.5*sys->dt / mvsq2e * part[i].fz / sys->mass;
//         part[i].rx += sys->dt*part[i].vx;
//         part[i].ry += sys->dt*part[i].vy;
//         part[i].rz += sys->dt*part[i].vz;
//     }
// }

// void verlet2(mdsys_t *sys,particles_t *part)
// {
//     /* second part: propagate velocities by another half step */

//     // add openMP parallel loop
//     #ifdef _OPENMP
//     #pragma omp parallel for
//     #endif
//     for (int i=0; i<sys->natoms; ++i) {
//         part[i].vx += 0.5*sys->dt / mvsq2e * part[i].fx / sys->mass;
//         part[i].vy += 0.5*sys->dt / mvsq2e * part[i].fy / sys->mass;
//         part[i].vz += 0.5*sys->dt / mvsq2e * part[i].fz / sys->mass;
//     }
// }

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
        vel[i].vx += 0.5*sys->dt / mvsq2e * forces[i].fx / sys->mass;
        vel[i].vy += 0.5*sys->dt / mvsq2e * forces[i].fy / sys->mass;
        vel[i].vz += 0.5*sys->dt / mvsq2e * forces[i].fz / sys->mass;
        // Update velocities
        // vel[i].vx += tmp * forces[i].fx/sys->mass;
        // vel[i].vy += tmp * forces[i].fy/sys->mass;
        // vel[i].vz += tmp * forces[i].fz/sys->mass;
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
        vel[i].vx += 0.5*sys->dt / mvsq2e * forces[i].fx / sys->mass;
        vel[i].vy += 0.5*sys->dt / mvsq2e * forces[i].fy / sys->mass;
        vel[i].vz += 0.5*sys->dt / mvsq2e * forces[i].fz / sys->mass;
        // vel[i].vx += tmp * forces[i].fx/sys->mass;
        // vel[i].vy += tmp * forces[i].fy/sys->mass;
        // vel[i].vz += tmp * forces[i].fz/sys->mass;
    }
}
