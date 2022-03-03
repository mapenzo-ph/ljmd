#include "utils.h"
#include "stdio.h"

#ifndef USE_MPI
/* append data to output. */
void output(mdsys_t *sys,particles_t *part, FILE *erg, FILE *traj)
{
    int i;

    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", part[i].rx, part[i].ry, part[i].rz);
    }
}
#else
extern int rank, size, myq;
extern MPI_Datatype MPI_PAR;
/* append data to output. 
        ACTUALLY THIS PART AS THE INPUT PART SHOULD BE DONE WITH MPI-IO OR HDF5 /!\ 
        WE DIDN'T HAVE THE TIME TO CHANGE THEM YET BUT THIS CAN BE DONE IN THE FUTURE!
*/
void output(mdsys_t *sys,particles_t *part, FILE *erg, FILE *traj, int *disp, int *count)
{
    int i;
    int cmode;
    MPI_Status status;
    if(rank==0){
        printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
        fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
        fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
        particles_t *buf=malloc(sizeof(particles_t)*sys->natoms);
        MPI_Gatherv(part,myq,MPI_PAR,buf,count,disp,MPI_PAR,0,MPI_COMM_WORLD);
        for (i=0; i<sys->natoms; ++i) {
                fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", buf[i].rx, buf[i].ry, buf[i].rz);
        }
    }else{
        MPI_Gatherv(part,myq,MPI_PAR,NULL,NULL,NULL,MPI_PAR,0,MPI_COMM_WORLD);
    }
}
#endif


