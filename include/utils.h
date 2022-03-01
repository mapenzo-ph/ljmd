#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
static const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
static const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

// force_compute.c
void force(mdsys_t *);
// verlet.c
void verlet1(mdsys_t *);
void verlet2(mdsys_t *);
// kinetic.c
void ekin(mdsys_t *);
// utils.c
double wallclock();
void azzero(double *, const int );
double pbc(double , const double );
int get_a_line(FILE *, char *);
// output.c
void output(mdsys_t *, FILE *, FILE *);

#endif
