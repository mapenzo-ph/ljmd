#include "utils.h"
#include <math.h>

/* compute forces */
void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    double epot = 0.0;
    int i,j;

    /* zero energy and forces */
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i,j,rx,ry,rz,r,ffac) reduction(+:epot)
    #endif
    for(i=0; i < (sys->natoms)-1; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);
                

                epot += 4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

            
                #ifdef _OPENMP
                #pragma omp atomic
                sys->fx[i] += rx / r * ffac;
                #pragma omp atomic
                sys->fy[i] += ry / r * ffac;
                #pragma omp atomic
                sys->fz[i] += rz / r * ffac;
                #pragma omp atomic
                sys->fx[j] -= rx / r * ffac;
                #pragma omp atomic
                sys->fy[j] -= ry / r * ffac;
                #pragma omp atomic
                sys->fz[j] -= rz / r * ffac;
                #else
                sys->fx[i] += rx / r * ffac;
                sys->fy[i] += ry / r * ffac;
                sys->fz[i] += rz / r * ffac;
                sys->fx[j] -= rx / r * ffac;
                sys->fy[j] -= ry / r * ffac;
                sys->fz[j] -= rz / r * ffac;
                #endif
            }
        }
    }

    sys->epot = epot;
}
