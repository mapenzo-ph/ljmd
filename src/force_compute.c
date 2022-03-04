#include "utils.h"
#include <math.h>

#ifndef USE_MPI
/* compute forces */
void force(mdsys_t *sys,coords_t * coord, for_t* forces)
{
    double ffac;
    double rx,ry,rz;
    double epot=0.0 ;
    int i,j;
// this change of paprameters done for avoiding expensive math like power, sqrt, division.
     double sigma1=sys->sigma;
     double sigma6=(sigma1)*(sigma1)*(sigma1)*(sigma1)*(sigma1)*(sigma1);
     double sigma12=sigma6*sigma6;
     double rsq,rinv,r6;
     double c12 = 4.0 * sys->epsilon * sigma6 * sigma6 ;
     double c6  = 4.0 * sys->epsilon * sigma6 ;
     double rcsq = sys->rcut * sys->rcut;

    /* zero energy and forces */
    azzero(forces, sys->natoms);


    #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i,j,rx,ry,rz,r6,ffac,rinv,rsq) reduction(+:epot)
    #endif
    for(i=0; i < (sys->natoms)-1; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* get distance between particle i and j */
            rx=pbc(coord[i].rx - coord[j].rx, 0.5*sys->box);
            ry=pbc(coord[i].ry - coord[j].ry, 0.5*sys->box);
            rz=pbc(coord[i].rz - coord[j].rz, 0.5*sys->box);
            rsq =rx*rx + ry*ry + rz*rz;


            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
              rinv=1.0/rsq;
              r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;

                epot += r6*(c12*r6-c6);


                #ifdef _OPENMP
                #pragma omp atomic
                forces[i].fx += rx * ffac;
                #pragma omp atomic
                forces[i].fy += ry  * ffac;
                #pragma omp atomic
                forces[i].fz += rz  * ffac;
                #pragma omp atomic
                forces[j].fx -= rx  * ffac;
                #pragma omp atomic
                forces[j].fy -= ry * ffac;
                #pragma omp atomic
                forces[j].fz -= rz * ffac;
                #else
                forces[i].fx += rx  * ffac;
                forces[i].fy += ry * ffac;
                forces[i].fz += rz  * ffac;
                forces[j].fx -= rx * ffac;
                forces[j].fy -= ry  * ffac;
                forces[j].fz -= rz * ffac;
                #endif
            }
        }
    }

    sys->epot = epot;
}
#else

/* compute forces */
void force(mdsys_t *sys, coords_t * coord, for_t* forces, for_t* forces2, MPI_Datatype MPI_FORCE,MPI_Datatype MPI_COORD, MPI_Op SUM_F)
{
    double r=0.0,ffac=0.0;
    double rx=0.0,ry=0.0,rz=0.0;
    int i,j,ii;

    /* zero energy and forces */
    double epot_temp=0.0;
    sys->epot=0.0;
    // azzero(forces2, sys->natoms);
    memset(forces2, 0.0, sys->natoms * sizeof(for_t));

    double sigma1=sys->sigma;
    double sigma6=(sigma1)*(sigma1)*(sigma1)*(sigma1)*(sigma1)*(sigma1);
    double sigma12=sigma6*sigma6;
    double rsq=0.0,rinv=0.0,r6=0.0;
    double c12 = 4.0 * sys->epsilon * sigma6 * sigma6 ;
    double c6  = 4.0 * sys->epsilon * sigma6 ;
    double rcsq = sys->rcut * sys->rcut;
    double bbb=0.5*sys->box;
    
    // can be avoid
    MPI_Bcast(coord,sys->natoms,MPI_COORD,0,MPI_COMM_WORLD);


    #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(ii,i,j,rx,ry,rz,r6,ffac,rinv,rsq) reduction(+:epot_temp)
    #endif
    for(i=0; i < (sys->natoms)-1; i+=sys->size) {
        ii=i+sys->rank;
        if(ii<=(sys->natoms-1)){
        for(j=ii+1; j < (sys->natoms); ++j) {

            /* get distance between particle i and j */
            rx=pbc(coord[ii].rx - coord[j].rx, bbb);
            ry=pbc(coord[ii].ry - coord[j].ry, bbb);
            rz=pbc(coord[ii].rz - coord[j].rz, bbb);
            rsq =rx*rx + ry*ry + rz*rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                rinv=1.0/rsq;
                r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                epot_temp += r6*c12*r6 - r6*c6;


                #ifdef _OPENMP
                #pragma omp atomic
                forces2[ii].fx += rx * ffac;
                #pragma omp atomic
                forces2[ii].fy += ry  * ffac;
                #pragma omp atomic
                forces2[ii].fz += rz  * ffac;
                #pragma omp atomic
                forces2[j].fx -= rx  * ffac;
                #pragma omp atomic
                forces2[j].fy -= ry * ffac;
                #pragma omp atomic
                forces2[j].fz -= rz * ffac;
                #else
                forces2[ii].fx += rx  * ffac;
                forces2[ii].fy += ry * ffac;
                forces2[ii].fz += rz  * ffac;
                forces2[j].fx -= rx * ffac;
                forces2[j].fy -= ry  * ffac;
                forces2[j].fz -= rz * ffac;
                #endif
            }
        }
        }
    }

    MPI_Reduce(forces2,forces,sys->natoms,MPI_FORCE,SUM_F,0,MPI_COMM_WORLD);
    MPI_Reduce(&epot_temp,&sys->epot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

}

#endif