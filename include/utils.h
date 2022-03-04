#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#ifdef USE_MPI
#include "mpi.h"
#include "string.h"

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
#ifndef USE_MPI
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
#else
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
    int size;
    int rank;
};
#endif

// struct _particles{
//     double rx,ry,rz;
//     double vx,vy,vz;
//     double fx,fy,fz;
// };

struct _coordinates{
    double rx,ry,rz;
};

struct _velocities{
    double vx,vy,vz;
};
struct _forces{
    double fx,fy,fz;
};

// typedef struct _particles particles_t;
typedef struct _mdsys mdsys_t;
typedef struct _coordinates coords_t;
typedef struct _velocities vel_t;
typedef struct _forces for_t;


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
#ifdef USE_MPI
// mpi_utils
MPI_Datatype mdsys_mpitype();
MPI_Datatype coordinates_mpitype();
MPI_Datatype velocities_mpitype();
MPI_Datatype forces_mpitype();
MPI_Op mpi_operation();
void sum_struct_force(void*, void* , int *, MPI_Datatype *);
#endif

// changed functions
void ekin(mdsys_t *, vel_t *);
void verlet1(mdsys_t *, coords_t *, vel_t *, for_t*);
void verlet2(mdsys_t *, vel_t *, for_t*);
#ifdef USE_MPI
void force(mdsys_t *, coords_t *, for_t*,for_t* , MPI_Datatype,MPI_Datatype, MPI_Op);
#else
void force(mdsys_t *, coords_t *, for_t*);
#endif
void azzero(for_t *, const int );
void output(mdsys_t *,coords_t *, FILE *, FILE *);

// void ekin(mdsys_t *, particles_t *);
// void verlet1(mdsys_t *, particles_t *);
// void verlet2(mdsys_t *, particles_t *);
// void force(mdsys_t *, particles_t *);
// void azzero(particles_t *, const int );

#endif
