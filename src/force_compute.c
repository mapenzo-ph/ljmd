#include "utils.h"
#include <math.h>

#ifndef USE_MPI
/* compute forces */
void force(mdsys_t *sys, particles_t *part)
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
    azzero(part, sys->natoms);


    #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i,j,rx,ry,rz,r6,ffac,rinv,rsq) reduction(+:epot)
    #endif
    for(i=0; i < (sys->natoms)-1; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* get distance between particle i and j */
            rx=pbc(part[i].rx - part[j].rx, 0.5*sys->box);
            ry=pbc(part[i].ry - part[j].ry, 0.5*sys->box);
            rz=pbc(part[i].rz - part[j].rz, 0.5*sys->box);
            rsq =rx*rx + ry*ry + rz*rz;


            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
              rinv=1.0/rsq;
              r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;

                epot += r6*(c12*r6-c6);


                #ifdef _OPENMP
                #pragma omp atomic
                part[i].fx += rx * ffac;
                #pragma omp atomic
                part[i].fy += ry  * ffac;
                #pragma omp atomic
                part[i].fz += rz  * ffac;
                #pragma omp atomic
                part[j].fx -= rx  * ffac;
                #pragma omp atomic
                part[j].fy -= ry * ffac;
                #pragma omp atomic
                part[j].fz -= rz * ffac;
                #else
                part[i].fx += rx  * ffac;
                part[i].fy += ry * ffac;
                part[i].fz += rz  * ffac;
                part[j].fx -= rx * ffac;
                part[j].fy -= ry  * ffac;
                part[j].fz -= rz * ffac;
                #endif
            }
        }
    }

    sys->epot = epot;
}
#else

extern int rank, size, myq;
extern MPI_Datatype MPI_PAR;
/* compute forces */
void force(mdsys_t *sys, particles_t *part)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(part, myq);


    //For now, later can be optimize
    int n_to_alloc=0;
    MPI_Allreduce(&myq,&n_to_alloc,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    particles_t *buf=calloc(sizeof(particles_t),n_to_alloc);
    particles_t *cbuf=calloc(sizeof(particles_t),n_to_alloc);

    int myq2;
    if(size>1)
    myq2=message(part,buf,myq);
    for(int l=0; l<size; ++l){
       for(i=0; i<n_to_alloc; ++i){ //nt_to_alloc and not myq because, otherwise some processors could not enter and a deadlock happens.
            for(j=0; j < myq2; ++j) {
                /* particles have no interactions with themselves.
                   this condition was substituted with r!=0, maybe we will change with r < tol=1e-12 to prevent machine error */
                // if (i==j && rank==l) continue;  

                /* get distance between particle i and j */
                rx=pbc(part[i].rx - buf[j].rx, 0.5*sys->box);
                ry=pbc(part[i].ry - buf[j].ry, 0.5*sys->box);
                rz=pbc(part[i].rz - buf[j].rz, 0.5*sys->box);
                r = sqrt(rx*rx + ry*ry + rz*rz);
                /* compute force and energy if within cutoff */
                if (r!=0 && r < sys->rcut) {
                    ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                            +6*pow(sys->sigma/r,6.0)/r);

                    sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                                -pow(sys->sigma/r,6.0));

                    part[i].fx += rx/r*ffac;
                    part[i].fy += ry/r*ffac;
                    part[i].fz += rz/r*ffac;
                }
            }
        }
        memcpy(cbuf,buf,sizeof(particles_t)*n_to_alloc);
        if(size>1)
        myq2=message(cbuf,buf,myq2);
    }

    MPI_Allreduce(MPI_IN_PLACE,&sys->epot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // TO DEBUG
    // for(int i=0; i<myq; ++i){
    //     printf("rank %d %lf %lf %lf \n",rank,part[i].fx,part[i].fy,part[i].fz);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("All done rank %d !\n",rank);
    // MPI_Abort(MPI_COMM_WORLD,0);
    free(buf);
    free(cbuf);

}


int message(particles_t *vec,particles_t *buf,int myq2){

    int q_send_recv[2]={myq2,0};

    if(rank==0){
        MPI_Send(&q_send_recv[0],1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
        MPI_Recv(&q_send_recv[1],1,MPI_INT,size-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else if(rank==size-1){
        MPI_Recv(&q_send_recv[1],1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&q_send_recv[0],1,MPI_INT,0,0,MPI_COMM_WORLD);
    }else{
        MPI_Recv(&q_send_recv[1],1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&q_send_recv[0],1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
    }

    // for(int i=0; i<myq2; ++i){
    //     printf("Before rank %d %lf %lf %lf \n",rank,vec[i].rx,vec[i].ry,vec[i].rz);
    // }
    // printf("After rank %d qsend %d qrecv %d \n",rank,q_send_recv[0],q_send_recv[1]);
    //     printf("Enter Here with rank %d !\n",rank);
    //     MPI_Abort(MPI_COMM_WORLD,0);


    if(rank==0){
        MPI_Send(vec,q_send_recv[0],MPI_PAR,rank+1,0,MPI_COMM_WORLD);
        MPI_Recv(buf,q_send_recv[1],MPI_PAR,size-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else if(rank==size-1){
        MPI_Recv(buf,q_send_recv[1],MPI_PAR,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(vec,q_send_recv[0],MPI_PAR,0,0,MPI_COMM_WORLD);
    }else{
        MPI_Recv(buf,q_send_recv[1],MPI_PAR,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(vec,q_send_recv[0],MPI_PAR,rank+1,0,MPI_COMM_WORLD);
    }

    // for(int i=0; i<q_send_recv[1]; ++i){
    //     printf("after rank %d %lf %lf %lf \n",rank,buf[i].rx,buf[i].ry,buf[i].rz);
    // }
    // MPI_Abort(MPI_COMM_WORLD,0);
    // MPI_Barrier(MPI_COMM_WORLD);


return q_send_recv[1];
}

#endif