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
void azzero(for_t *d, const int n)
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
MPI_Datatype coordinates_mpitype(){
    MPI_Datatype COORD_TYPE;
    MPI_Aint disp[3];
    MPI_Aint base_addr;
    int length[3]={1,1,1};
    coords_t tmp;
    MPI_Get_address(&tmp,&base_addr);
    MPI_Get_address(&tmp.rx,&disp[0]);
    MPI_Get_address(&tmp.ry,&disp[1]);
    MPI_Get_address(&tmp.rz,&disp[2]);
    disp[0] = MPI_Aint_diff(disp[0], base_addr);
    disp[1] = MPI_Aint_diff(disp[1], base_addr);
    disp[2] = MPI_Aint_diff(disp[2], base_addr);
    MPI_Datatype type[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Type_create_struct(3, length, disp, type, &COORD_TYPE);
    MPI_Type_commit(&COORD_TYPE);
    return COORD_TYPE;
}

MPI_Datatype velocities_mpitype(){
    MPI_Datatype VEL_TYPE;
    MPI_Aint disp[3];
    MPI_Aint base_addr;
    int length[3]={1,1,1};
    vel_t tmp;
    MPI_Get_address(&tmp,&base_addr);
    MPI_Get_address(&tmp.vx,&disp[0]);
    MPI_Get_address(&tmp.vy,&disp[1]);
    MPI_Get_address(&tmp.vz,&disp[2]);
    disp[0] = MPI_Aint_diff(disp[0], base_addr);
    disp[1] = MPI_Aint_diff(disp[1], base_addr);
    disp[2] = MPI_Aint_diff(disp[2], base_addr);
    MPI_Datatype type[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Type_create_struct(3, length, disp, type, &VEL_TYPE);
    MPI_Type_commit(&VEL_TYPE);
    return VEL_TYPE;
}

MPI_Datatype forces_mpitype(){
    MPI_Datatype FORCE_TYPE;
    MPI_Aint disp[3];
    MPI_Aint base_addr;
    int length[3]={1,1,1};
    for_t tmp;
    MPI_Get_address(&tmp,&base_addr);
    MPI_Get_address(&tmp.fx,&disp[0]);
    MPI_Get_address(&tmp.fy,&disp[1]);
    MPI_Get_address(&tmp.fz,&disp[2]);
    disp[0] = MPI_Aint_diff(disp[0], base_addr);
    disp[1] = MPI_Aint_diff(disp[1], base_addr);
    disp[2] = MPI_Aint_diff(disp[2], base_addr);
    MPI_Datatype type[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Type_create_struct(3, length, disp, type, &FORCE_TYPE);
    MPI_Type_commit(&FORCE_TYPE);
    return FORCE_TYPE;
}

// Ibuf -> input buffer
// Obuf -> output buffer
// MPI_Datatype
void sum_struct_force(void* Ibuf, void* Obuf, int *len, MPI_Datatype *MPI_DATA){
    for_t* input = (for_t*)Ibuf;
    for_t* output = (for_t*)Obuf;

    for(int i = 0; i < *len; i++)
    {
        output[i].fx += input[i].fx;
        output[i].fy += input[i].fy;
        output[i].fz += input[i].fz;
    }
}

MPI_Op mpi_operation(){
    //Creation operation for force compute
    MPI_Op operation;
    MPI_Op_create(&sum_struct_force, 1, &operation);
    return operation;
}
#endif
