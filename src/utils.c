#include "utils.h"
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* helper function: get current time in seconds since epoch */
double wallclock()
{
        struct timeval t;
        gettimeofday(&t,0);
        return ((double) t.tv_sec) + 1.0e-6*((double) t.tv_usec);
}

/* helper function: zero out an array */
void azzero(particles_t *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i].fx=d[i].fy=d[i].fz=0.0;
    }
}

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;
        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}

#ifdef USE_MPI
/* helper function: for mpi implementation
                    generates mpi datastructure for _mdsys structure  */
// This could be simplify 
MPI_Datatype mdsys_mpitype(){
    MPI_Datatype MDSYS_TYPE;
    MPI_Aint disp[12];
    MPI_Aint base_addr;
    int length[12]={1,1,1,1,1,1,1,1,1,1,1,1};
    mdsys_t tmp;
    MPI_Get_address(&tmp,&base_addr);
    MPI_Get_address(&tmp.dt,&disp[0]);
    MPI_Get_address(&tmp.mass,&disp[1]);
    MPI_Get_address(&tmp.epsilon,&disp[2]);
    MPI_Get_address(&tmp.sigma,&disp[3]);
    MPI_Get_address(&tmp.box,&disp[4]);
    MPI_Get_address(&tmp.rcut,&disp[5]);
    MPI_Get_address(&tmp.ekin,&disp[6]);
    MPI_Get_address(&tmp.epot,&disp[7]);
    MPI_Get_address(&tmp.temp,&disp[8]);
    MPI_Get_address(&tmp.natoms,&disp[9]);
    MPI_Get_address(&tmp.nfi,&disp[10]);
    MPI_Get_address(&tmp.nsteps,&disp[11]);

    disp[0] = MPI_Aint_diff(disp[0], base_addr);
    disp[1] = MPI_Aint_diff(disp[1], base_addr);
    disp[2] = MPI_Aint_diff(disp[2], base_addr);
    disp[3] = MPI_Aint_diff(disp[3], base_addr);
    disp[4] = MPI_Aint_diff(disp[4], base_addr);
    disp[5] = MPI_Aint_diff(disp[5], base_addr);
    disp[6] = MPI_Aint_diff(disp[6], base_addr);
    disp[7] = MPI_Aint_diff(disp[7], base_addr);
    disp[8] = MPI_Aint_diff(disp[8], base_addr);
    disp[9] = MPI_Aint_diff(disp[9], base_addr);
    disp[10] = MPI_Aint_diff(disp[10], base_addr);
    disp[11] = MPI_Aint_diff(disp[11], base_addr);

    MPI_Datatype type[12]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                        MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                        MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                        MPI_INT,MPI_INT,MPI_INT};
    
    MPI_Type_create_struct(12, length, disp, type, &MDSYS_TYPE);
    MPI_Type_commit(&MDSYS_TYPE);
    return MDSYS_TYPE;
}

MPI_Datatype particle_mpitype(){
    MPI_Datatype PAR_TYPE;
    MPI_Aint disp[9];
    MPI_Aint base_addr;
    int length[9]={1,1,1,1,1,1,1,1,1};
    particles_t tmp;
    MPI_Get_address(&tmp,&base_addr);
    MPI_Get_address(&tmp.rx,&disp[0]);
    MPI_Get_address(&tmp.ry,&disp[1]);
    MPI_Get_address(&tmp.rz,&disp[2]);
    MPI_Get_address(&tmp.vx,&disp[3]);
    MPI_Get_address(&tmp.vy,&disp[4]);
    MPI_Get_address(&tmp.vz,&disp[5]);
    MPI_Get_address(&tmp.fx,&disp[6]);
    MPI_Get_address(&tmp.fy,&disp[7]);
    MPI_Get_address(&tmp.fz,&disp[8]);

    disp[0] = MPI_Aint_diff(disp[0], base_addr);
    disp[1] = MPI_Aint_diff(disp[1], base_addr);
    disp[2] = MPI_Aint_diff(disp[2], base_addr);
    disp[3] = MPI_Aint_diff(disp[3], base_addr);
    disp[4] = MPI_Aint_diff(disp[4], base_addr);
    disp[5] = MPI_Aint_diff(disp[5], base_addr);
    disp[6] = MPI_Aint_diff(disp[6], base_addr);
    disp[7] = MPI_Aint_diff(disp[7], base_addr);
    disp[8] = MPI_Aint_diff(disp[8], base_addr);
    MPI_Datatype type[9]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                    MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                    MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Type_create_struct(9, length, disp, type, &PAR_TYPE);
    MPI_Type_commit(&PAR_TYPE);
    return PAR_TYPE;
}
#endif
