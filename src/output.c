#include "utils.h"
#include "stdio.h"

/* append data to output. */
void output(mdsys_t *sys,coords_t *coord, FILE *erg, FILE *traj)
{
    int i;
    #ifdef USE_MPI
    if(sys->rank==0){
        printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
        fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
        fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
        for (i=0; i<sys->natoms; ++i) {
                fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", coord[i].rx, coord[i].ry, coord[i].rz);
        }
    }
    #else
        printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
        fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
        fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
        for (i=0; i<sys->natoms; ++i) {
                fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", coord[i].rx, coord[i].ry, coord[i].rz);
        }
    #endif
}


