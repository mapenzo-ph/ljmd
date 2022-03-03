#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#ifdef USE_MPI
#include "mpi.h"
#endif

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
    double dt;
    double mass;
    double epsilon;
    double sigma;
    double box;
    double rcut;
    double ekin;
    double epot;
    double temp;
    int natoms;
    int nfi;
    int nsteps;
};

struct _particles{
    double rx,ry,rz;
    double vx,vy,vz;
    double fx,fy,fz;
};

typedef struct _mdsys mdsys_t;
typedef struct _particles particles_t;
// // force_compute.c
// void force(mdsys_t *, particles_t *);
// // verlet.c
// void verlet1(mdsys_t *);
// void verlet2(mdsys_t *);
// // kinetic.c
// void ekin(mdsys_t *);
// // utils.c
double wallclock();
// void azzero(double *, const int );
double pbc(double , const double );
int get_a_line(FILE *, char *);
// output.c
// MPI
// changed functions
void ekin(mdsys_t *, particles_t *);
void verlet1(mdsys_t *, particles_t *);
void verlet2(mdsys_t *, particles_t *);
void force(mdsys_t *, particles_t *);
void azzero(particles_t *, const int );

#ifndef USE_MPI
void output(mdsys_t *,particles_t *, FILE *, FILE *);
#else
// global variable
int rank, size;
int myq;
MPI_Datatype MPI_PAR;
void output(mdsys_t *,particles_t *, FILE *, FILE *, int *, int *);
// message
int message(particles_t *,particles_t *, int );
MPI_Datatype mdsys_mpitype();
MPI_Datatype particle_mpitype();
#endif

#endif
