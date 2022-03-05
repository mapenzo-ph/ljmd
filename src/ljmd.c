/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "utils.h"


/* main */
int main(int argc, char **argv)
{
    mdsys_t sys;

    #ifdef USE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&sys.rank);
    MPI_Comm_size(MPI_COMM_WORLD,&sys.size);
    // --------------------------------------------------------------------------
    // Declaration of some structures 
    coords_t* coord=NULL;
    vel_t* vel=NULL;
    for_t* forces=NULL;
    for_t* forces2=NULL;
    // --------------------------------------------------------------------------
    // Creation common datatype
    MPI_Datatype MPI_COORD=coordinates_mpitype();
    // MPI_VEL: Not use in this parallelization strategy but maybe useful for another strategy
    // MPI_Datatype MPI_VEL=velocities_mpitype();  
    MPI_Datatype MPI_FORCE=coordinates_mpitype();
    MPI_Op MPI_SUM_F=mpi_operation();
    // --------------------------------------------------------------------------
    #else
    coords_t* coord=NULL;
    vel_t* vel=NULL;
    for_t* forces=NULL;
    #endif


    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *traj,*erg,*fp,*gh;  
    double t_start;

    #ifdef USE_MPI
        if(sys.rank==0){
            printf("LJMD version %3.1f\n", LJMD_VERSION);
            t_start = MPI_Wtime();
        }
    #else
        printf("LJMD version %3.1f\n", LJMD_VERSION);
        t_start = wallclock();
    #endif

    gh=fopen(argv[1],"r");

    /* read input file */
    if(get_a_line(gh,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(gh,line)) return 1;
    sys.mass=atof(line);
    if(get_a_line(gh,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(gh,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(gh,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(gh,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(gh,restfile)) return 1;
    if(get_a_line(gh,trajfile)) return 1;
    if(get_a_line(gh,ergfile)) return 1;
    if(get_a_line(gh,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(gh,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(gh,line)) return 1;
    nprint=atoi(line);

    #ifdef USE_MPI
        coord=calloc(sizeof(coords_t),sys.natoms);
        forces2=(for_t*)malloc(sizeof(for_t)*sys.natoms);
        if(sys.rank==0){
            vel=calloc(sizeof(vel_t),sys.natoms);
            forces=(for_t*)malloc(sizeof(for_t)*sys.natoms);
            /* read restart */
            fp=fopen(restfile,"r");
            if(fp) {
                for (i=0; i<sys.natoms; ++i) {
                    fscanf(fp,"%lf%lf%lf",&coord[i].rx, &coord[i].ry, &coord[i].rz);
                }
                for (i=0; i<sys.natoms; ++i) {
                    fscanf(fp,"%lf%lf%lf",&vel[i].vx, &vel[i].vy, &vel[i].vz);
                }
                fclose(fp);
                azzero(forces, sys.natoms);
            } else {
                perror("cannot read restart file");
                return 3;
            }
        }
    #else
        vel=calloc(sizeof(vel_t),sys.natoms);
        forces=(for_t*)malloc(sizeof(for_t)*sys.natoms);
        coord=calloc(sizeof(coords_t),sys.natoms);
        /* read restart */
        fp=fopen(restfile,"r");
        if(fp) {
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",&coord[i].rx, &coord[i].ry, &coord[i].rz);
            }
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",&vel[i].vx, &vel[i].vy, &vel[i].vz);
            }
            fclose(fp);
            azzero(forces, sys.natoms);
        } else {
            perror("cannot read restart file");
            return 3;
        }
    #endif

    /* initialize forces and energies.*/
    #ifdef USE_MPI
    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys,coord,forces,forces2,MPI_FORCE,MPI_COORD,MPI_SUM_F);
    if(sys.rank==0)
        ekin(&sys,vel);
    #else 
    sys.nfi=0;
    force(&sys,coord,forces);
    ekin(&sys,vel);
    #endif

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    #ifdef USE_MPI
        double t_end=MPI_Wtime()-t_start;
        MPI_Allreduce(MPI_IN_PLACE,&t_end,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        if(sys.rank==0){
            printf("Startup time: %10.3fs\n", t_end);
            printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
            printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        }
        output(&sys,coord, erg, traj);
        t_start = MPI_Wtime();
    #else
        printf("Startup time: %10.3fs\n", wallclock()-t_start);
        printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys,coord, erg, traj);
        /* reset timer */
        t_start = wallclock();
    #endif


    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys,coord, erg, traj);

        #ifdef USE_MPI
            if(sys.rank==0)
                verlet1(&sys,coord,vel,forces);
            force(&sys,coord,forces,forces2,MPI_FORCE,MPI_COORD,MPI_SUM_F);
            if(sys.rank==0){
                verlet2(&sys,vel,forces); 
                ekin(&sys,vel);
            }
        #else
            /* propagate system and recompute energies */
            verlet1(&sys,coord,vel,forces);
            force(&sys,coord,forces);
            verlet2(&sys,vel,forces); 
            ekin(&sys,vel);
        #endif

    }
    /**************************************************/

    /* clean up: close files, free memory */
    #ifdef USE_MPI
    t_end=MPI_Wtime()-t_start;
    MPI_Allreduce(MPI_IN_PLACE,&t_end,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(sys.rank==0){
        printf("Simulation Done. Run time: %10.3fs\n", t_end);
    }
    #else
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
    #endif
    fclose(erg);
    fclose(traj);

    #ifdef USE_MPI
    free(coord);
    free(vel);
    free(forces);
    free(forces2);
    MPI_Finalize();
    #else
    free(coord);
    free(vel);
    free(forces);
    #endif

    return 0;
}


