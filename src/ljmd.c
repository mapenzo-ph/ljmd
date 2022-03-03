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

#ifndef USE_MPI

/* main */
int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);

    // exit(0);

    /* allocate memory */

    particles_t *part=calloc(sizeof(particles_t),sys.natoms);

    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",&part[i].rx, &part[i].ry, &part[i].rz);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",&part[i].vx, &part[i].vy, &part[i].vz);
        }
        fclose(fp);
        azzero(part, sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys,part);
    ekin(&sys,part);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Startup time: %10.3fs\n", wallclock()-t_start);
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys,part, erg, traj);

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys,part, erg, traj);

        /* propagate system and recompute energies */
        verlet1(&sys,part);
        force(&sys,part);

        verlet2(&sys,part);
        ekin(&sys,part);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
    fclose(erg);
    fclose(traj);

    free(part);

    return 0;
}

#else
// global variable
int rank, size, myq;
MPI_Datatype MPI_PAR;

/* main */
int main(int argc, char **argv)
{

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    mdsys_t sys;
    particles_t *part=NULL;
    MPI_PAR=particle_mpitype();
    int nprint;
    int count[size];
    int disp[size];
    char trajfile[BLEN], ergfile[BLEN];
    FILE *traj,*erg;
    if(rank==0){
        int i;
        char restfile[BLEN], line[BLEN];
        FILE *fp;  
        // printf("LJMD version %3.1f\n", LJMD_VERSION);
        printf("LJMD version %3.1f\n", 1.0);

        /* read input file */
        if(get_a_line(stdin,line)) return 1;
        sys.natoms=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.mass=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.epsilon=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.sigma=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.rcut=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.box=atof(line);
        if(get_a_line(stdin,restfile)) return 1;
        if(get_a_line(stdin,trajfile)) return 1;
        if(get_a_line(stdin,ergfile)) return 1;
        if(get_a_line(stdin,line)) return 1;
        sys.nsteps=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.dt=atof(line);
        if(get_a_line(stdin,line)) return 1;
        nprint=atoi(line);

        particles_t *tmp=calloc(sizeof(particles_t),sys.natoms);

        /* read restart */
        fp=fopen(restfile,"r");
        if(fp) {
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",&tmp[i].rx, &tmp[i].ry, &tmp[i].rz);
            }
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",&tmp[i].vx, &tmp[i].vy, &tmp[i].vz);
            }
            fclose(fp);
            azzero(tmp, sys.natoms);
        } else {
            perror("cannot read restart file");
            return 3;
        }

        //  set up particles for each cores
        myq=sys.natoms/size;
        int mod=sys.natoms%size;
        for(int i=0; i<size; ++i){
            if(i<mod){
                count[i]=myq+1;
                disp[i]=i*count[i];
            }else{
                count[i]=myq;
                disp[i]=i*count[i]+mod;
            }
        }

        erg=fopen(ergfile,"w");
        traj=fopen(trajfile,"w");

        // TO DEBUG
        // if(size==4){
        //     printf("count %d %d %d %d\n",count[0],count[1],count[2],count[3]);
        //     printf("disp %d %d %d %d\n",disp[0],disp[1],disp[2],disp[3]);
        // }

        MPI_Scatter(count,1,MPI_INT,&myq,1,MPI_INT,0,MPI_COMM_WORLD);
        part=malloc(sizeof(particles_t)*myq);
        MPI_Scatterv(tmp,count,disp,MPI_PAR,part,myq,MPI_PAR,0,MPI_COMM_WORLD);
        free(tmp);
    }else{
        MPI_Scatter(NULL,1,MPI_INT,&myq,1,MPI_INT,0,MPI_COMM_WORLD);
        part=malloc(sizeof(particles_t)*myq);
        MPI_Scatterv(NULL,NULL,NULL,MPI_PAR,part,myq,MPI_PAR,0,MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Datatype MPI_SYS=mdsys_mpitype();
    MPI_Bcast(&sys,1,MPI_SYS,0,MPI_COMM_WORLD);

    //Broadcast file name- try failed for now
    MPI_Bcast(&nprint,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(disp,size,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(count,size,MPI_INT,0,MPI_COMM_WORLD);
    // MPI_FILE* fh,gh;
    // MPI_Bcast(ergfile,BLEN,MPI_CHAR,0,MPI_COMM_WORLD);
    // MPI_Bcast(trajfile,BLEN,MPI_CHAR,0,MPI_COMM_WORLD);

    // TO DEBUG
    // for(int i=0; i<myq; ++i){
    //     // printf("rank %d %lf %lf %lf \n",rank,part[i].rx,part[i].ry,part[i].rz);
    //     // printf("rank %d %lf %lf %lf \n",rank,part[i].vx,part[i].vy,part[i].vz);        
    // }

    //TO DEBUG
    // printf("rank %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d \n",rank,
    //     sys.dt,sys.mass,sys.epsilon,sys.sigma,sys.box,sys.rcut,sys.ekin,sys.epot,sys.temp,
    //     sys.natoms,sys.nfi,sys.nsteps);



    // /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys,part);
    ekin(&sys,part);

    if(rank==0){
    // printf("Startup time: %10.3fs\n", wallclock()-t_start);
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    }
    output(&sys,part,erg,traj,disp,count);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
                output(&sys,part,erg,traj,disp,count);
        /* propagate system and recompute energies */

        verlet1(&sys,part);
        force(&sys,part);
        verlet2(&sys,part);
        ekin(&sys,part);

    }
    /**************************************************/

    // /* clean up: close files, free memory */

    if(rank==0){    
        fclose(erg);
        fclose(traj);
    }
    free(part);

    MPI_Finalize();


    return 0;
}


#endif

