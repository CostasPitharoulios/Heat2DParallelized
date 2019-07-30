/****************************************************************************
 * FILE: mpi_heat2D.c
 * DESCRIPTIONS:  
 *   HEAT2D Example - Parallelized C Version
 *   This example is based on a simplified two-dimensional heat 
 *   equation domain decomposition.  The initial temperature is computed to be 
 *   high in the middle of the domain and zero at the boundaries.  The 
 *   boundaries are held at zero throughout the simulation.  During the 
 *   time-stepping, an array containing two domains is used; these domains 
 *   alternate between old data and new data.
 *
 *   In this parallelized version, the grid is decomposed by the master
 *   process and then distributed by rows to the worker processes.  At each 
 *   time step, worker processes must exchange border data with neighbors, 
 *   because a grid point's current temperature depends upon it's previous
 *   time step value plus the values of the neighboring grid points.  Upon
 *   completion of all time steps, the worker processes return their results
 *   to the master process.
 *
 *   Two data files are produced: an initial data set and a final data set.
 * AUTHOR: Blaise Barney - adapted from D. Turner's serial C version. Converted
 *   to MPI: George L. Gusciora (1/95)
 * LAST REVISED: 04/02/05
 ****************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NXPROB      12                 /* x dimension of problem grid */
#define NYPROB      8                 /* y dimension of problem grid */
#define STEPS       1 /*100*/            /* number of time steps */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};

void inidat(), prtdat(), update(), myprint(), DUMMYDUMDUM();
int malloc2darr(),free2darr(),isPrime();

int main (int argc, char *argv[]){
    float u[2][NXPROB][NYPROB],        /* array for grid TODO: mhpws na to exei mono o master? den xreiazetai na desmeutei se olous.. oi uloipoi exoyn to local. (auto mporei na ginei vazontas to static mesa se if, isws) */
          /* Episis den xreiazomaste u[2] efoson exoume local[2] */
          **local[2];                  /* array for local part of the grid */
    int	taskid,                     /* this task's unique id */
        numworkers,                 /* number of worker processes */
        numtasks,                   /* number of tasks */
        /*averow,offset,*/offsetX, offsetY,/*extra,*/   /* for sending rows of data */
        dest, source,               /* to - from for message send-receive */
        left,right,up,down,        /* neighbor tasks */
        msgtype,                    /* for message types */
        rc,start,end,               /* misc */
        xdim, ydim,                 /* dimensions of grid partition (e.x. 4x4) */
        rows, columns,             /* number of rows/columns of each block (e.x. 20x12) */
        i,j,x,y,ix,iy,iz,it;              /* loop variables */
    MPI_Status status;


    /* First, find out my taskid and how many tasks are running */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers = numtasks;

    if (taskid == MASTER) {
        /************************* Master code *******************************/

        if ((isPrime(numworkers))){
            printf("ERROR: the number of workers is prime (%d).\n",numworkers);
            MPI_Abort(MPI_COMM_WORLD, 22);
            exit(22);
        }
        printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

        /* If the number of cells is not divisible by numworkers, abort */
        if ((NXPROB*NYPROB)%numworkers){
            printf("ERROR: number of cells is not divisible by the number of workers\n");
            MPI_Abort(MPI_COMM_WORLD, 22);
            exit(22);
        }

        /* Initialize grid */
        printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
        printf("Initializing grid and writing initial.dat file...\n");
        DUMMYDUMDUM(NXPROB, NYPROB, u); /* TODO TODO TODO */
        prtdat(NXPROB, NYPROB, u, "initial.dat");
        for (ix=0; ix<NXPROB; ix++){
            for (j=0; j<NYPROB; j++)
                printf("%6.1f ", u[0][ix][j]);
            printf("\n\n");
        }
        //myprint(NXPROB, NYPROB, u[0]);

        /* Find the dimentions of the partitioned grid (e.x. 4 x 4) */
        /* xdim,ydim are guarented to be found, since we have checked that
         * numworkers is not prime. */
        for (x=sqrt(numworkers) + 1; x>=1; x--){
            if (numworkers % x == 0){
                xdim = x;
                ydim = numworkers/x;
                break;
            }
        }
       
        /* Swap them if neccessary, in order to make the blocks more square-like */ 
        if (NYPROB > NXPROB && ydim < xdim){
            int a = xdim;
            xdim = ydim;
            ydim = a;
        }

        printf("The grid will part into a %d x %d block grid.\n",xdim,ydim);

        /* Compute the length and height of each block */
        rows = NXPROB / xdim;
        columns = NYPROB / ydim;
        printf("Each block is %d x %d.\n",rows,columns);

        /* Distribute work to workers.*/ 
        //offsetX = 0;
        //offsetY = 0;
        for (i=1; i<numworkers; i++){
         ///rows = (i <= extra) ? averow+1 : averow; 

            /* Compute the coordinates of the up left corner of the block */
            //offsetX = ((i-1)%xdim)*columns; /* TODO isws na htan pio oikonomiko na ekmetaleutoume oti eimaste se for loop opws eipe o kwstas */
            //offsetY = ((i-1)/xdim)*rows;

            /* Find the neighbours of this block */
            if (i < ydim) // if this is the first row
                up = MPI_PROC_NULL;
            else
                up = i - ydim;

            if (i >= ((xdim-1) * ydim)) //if this is the last row
               down = MPI_PROC_NULL;
            else
               down = i + ydim;

            if (i%ydim == 0)	// if this is the first column
                left = MPI_PROC_NULL;
            else
                left = i-1;

            if (i%ydim == 3)	//if this is the last column
                right = MPI_PROC_NULL;
            else
                right = i+1;



            /*  Now send startup information to each worker  */
            /*TODO isws kai auta prepei na ginoun ISend */
            dest = i;
            //MPI_Send(&offsetX, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            //MPI_Send(&offsetY, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&columns, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&down, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            //MPI_Send(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);
            //printf("Sent to task %d: rows= %d offset= %d ",dest,rows,offset);
            //printf("left= %d right= %d\n",left,right);
            //offsetX = offsetX + columns;  //PEIRAKSA KAI AUTA NA EINAI ETOIMA
            //offsetY = offsetY + rows;	//PEIRAKSA KAI AUTA NA EINAI ETOIMA
        }

        /* Master does its part of the work */
        //offsetX = 0;
        //offsetY = 0;
        left = MPI_PROC_NULL;
        right = 1;
        up = MPI_PROC_NULL;
        down = ydim;

    }else{
        /*************** workers code *****************/

        /* Initialize everything - including the borders - to zero */
        for (iz=0; iz<2; iz++)
            for (ix=0; ix<NXPROB; ix++) 
                for (iy=0; iy<NYPROB; iy++) 
                    u[iz][ix][iy] = 0.0;

        /* Receive my offset, rows, neighbors and grid partition from master */
        source = MASTER; msgtype = BEGIN;
        //MPI_Recv(&offsetX, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        //MPI_Recv(&offsetY, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&columns, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    }
    printf("LOG: Process %d: left:%d, right:%d, up:%d, down:%d\n",taskid,left,right,up,down);

    malloc2darr(&local[0], rows+2, columns+2);
    malloc2darr(&local[1], rows+2, columns+2);

    /* Initialize with 0's TODO mhpws den xreiazetai? */
    for (iz=0; iz<2; iz++)
        for (ix=0; ix<rows+2; ix++) 
            for (iy=0; iy<columns+2; iy++) 
                local[iz][ix][iy] = 0.0;



    /* Define the datatype of send buffer elements */
    int sendsizes[2]    = {NXPROB, NYPROB};    /* u size */
    int sendsubsizes[2] = {rows, columns};     /* local size without halo */
    int sendstarts[2]   = {0,0};

    MPI_Datatype type, sendsubarrtype;
    MPI_Type_create_subarray(2, sendsizes, sendsubsizes, sendstarts, MPI_ORDER_C, MPI_FLOAT, &type);
    MPI_Type_create_resized(type, 0, columns*sizeof(float), &sendsubarrtype); /* h columns */
    MPI_Type_commit(&sendsubarrtype);

    /* Define the datatype of receive buffer elements */
    int recvsizes[2]    = {rows+2, columns+2};         /* local array size */
    int recvsubsizes[2] = {rows, columns};          /* local size without halo */
    int recvstarts[2]   = {0,0};

    MPI_Datatype recvsubarrtype;
    MPI_Type_create_subarray(2, recvsizes, recvsubsizes, recvstarts, MPI_ORDER_C, MPI_FLOAT, &recvsubarrtype);
    MPI_Type_commit(&recvsubarrtype);


    float *globalptr=NULL;
    if (taskid == MASTER) globalptr = &(u[0][0][0]);

    /* Scatter array to all processes */
    int *sendcounts = (int*)malloc(sizeof(int)*xdim*ydim);
    int *displs = (int*)malloc(sizeof(int)*xdim*ydim);
    
    if (taskid == MASTER){
        /* Every process has one piece */
        for (i=0; i<xdim*ydim; i++) sendcounts[i]=1; 

        /* Determine the starting point of every task's data */
        int disp = 0;
        for (i=0; i<xdim; i++){
            for (j=0; j<ydim; j++){
                //printf("displs[%d]=%d\n",i*ydim+j,disp);
                displs[i*ydim+j] = disp; /* h' mhpws xdim */
                disp +=1;
            }
            disp += (rows-1)*ydim; /* h' rows, ydim klp */
        }
    }

    if (taskid == 0){
        printf("displs=[ ");
        for (i=0; i<ydim*xdim; i++){
                printf("%d ",displs[i]);
        }
        printf("]\n");
    }


    MPI_Scatterv(globalptr, sendcounts, displs, sendsubarrtype, &(local[0][1][1]), columns*rows, recvsubarrtype, MASTER, MPI_COMM_WORLD);


    for ( i=0; i<numtasks; i++){
        if (taskid == i){
            printf("=========== To kommati tou %d =========\n",i);
            for (ix=0; ix<rows+2; ix++){
                for (j=0; j<columns+2; j++)
                    printf("%6.1f ", local[0][ix][j]);
                printf("\n\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

#if 0 
    if ( taskid!=MASTER){
    /* workers code */

      /* Determine border elements.  Need to consider first and last columns. */
      /* Obviously, row 0 can't exchange with row 0-1.  Likewise, the last */
      /* row can't exchange with last+1.  */
      if (offset==0) 
         start=1;
      else 
         start=offset;
      if ((offset+rows)==NXPROB) 
         end=start+rows-2;
      else 
         end = start+rows-1;

      /* Begin doing STEPS iterations.  Must communicate border rows with */
      /* neighbors.  If I have the first or last grid row, then I only need */
      /*  to  communicate with one neighbor  */
      printf("Task %d received work. Beginning time steps...\n",taskid);
      iz = 0;
      for (it = 1; it <= STEPS; it++)
      {
         if (left != NONE)
         {
            MPI_Send(&u[iz][offset][0], NYPROB, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
            source = left;
            msgtype = LTAG;
            MPI_Recv(&u[iz][offset-1][0], NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
         }
         if (right != NONE)
         {
            MPI_Send(&u[iz][offset+rows-1][0], NYPROB, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
            source = right;
            msgtype = RTAG;
            MPI_Recv(&u[iz][offset+rows][0], NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);
         }
         /* Now call update to update the value of grid points */
         update(start,end,NYPROB,&u[iz][0][0],&u[1-iz][0][0]);
         iz = 1 - iz;
      }
      
      /* Finally, send my portion of final results back to master */
      MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
      MPI_Send(&u[iz][offset][0], rows*NYPROB, MPI_FLOAT, MASTER, DONE, 
               MPI_COMM_WORLD);
   }
#endif

    /////////////////////
    /* Vazw ka8e diergasia na alla3ei ton local, gia testing */
    for (i=0; i<rows+2; i++){
        for (j=0; j<columns+2; j++){
            local[0][i][j] = taskid;
        }
    }
    /////////////////////

    /* Gather it all back */
    MPI_Gatherv(&(local[0][1][1]), columns*rows,  MPI_FLOAT/*recvsubarrtype*/, globalptr, sendcounts, displs, sendsubarrtype, 0, MPI_COMM_WORLD);

    free2darr(&local[0]);
    free2darr(&local[1]);

    MPI_Type_free(&sendsubarrtype); /*TODO kai tous upoloipous */

    if (taskid==MASTER){
        printf("Processed grid:\n");
        for (ix=0; ix<NXPROB; ix++){
            for (j=0; j<NYPROB; j++)
                printf("%6.1f ", u[0][ix][j]);
            printf("\n\n");
        }
    }

    MPI_Finalize();
    return 0;
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, float *u1, float *u2)
{
   int ix, iy;
   for (ix = start; ix <= end; ix++) 
      for (iy = 1; iy <= ny-2; iy++) 
         *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                          parms.cx * (*(u1+(ix+1)*ny+iy) +
                          *(u1+(ix-1)*ny+iy) - 
                          2.0 * *(u1+ix*ny+iy)) +
                          parms.cy * (*(u1+ix*ny+iy+1) +
                         *(u1+ix*ny+iy-1) - 
                          2.0 * *(u1+ix*ny+iy));
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
int ix, iy;

for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
int ix, iy;
FILE *fp;

fp = fopen(fnam, "w");
for (iy = ny-1; iy >= 0; iy--) {
  for (ix = 0; ix <= nx-1; ix++) {
    fprintf(fp, "%6.1f", *(u1+ix*ny+iy));
    if (ix != nx-1) 
      fprintf(fp, " ");
    else
      fprintf(fp, "\n");
    }
  }
fclose(fp);
}

/* Checkis if a given integer is a prime number */
int isPrime(int n){
    int i;
    if (n==2)
        return 1;
    if (n%2==0)
        return 0;
    for (i=3;i*i<=n;i+=2)
        if (n%i==0) 
            return 0;
    return 1;
}

/* TODO delete this func before we paradwsoume */
void myprint(int nx, int ny, float *u1) {
    int ix, iy;
    for (iy = ny-1; iy >= 0; iy--) {
      for (ix = 0; ix <= nx-1; ix++) {
        printf("%6.1f", *(u1+ix*ny+iy));
        if (ix != nx-1) 
          printf(" ");
        else
          printf("\n");
        }
    }
}



int malloc2darr(float ***array, int n, int m) {

    /* allocate the n*m contiguous items */
    float *p = (float *)malloc(n*m*sizeof(float));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (float **)malloc(n*sizeof(float*));
    if (!(*array)) {
        free(p);
        return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n; i++)
        (*array)[i] = &(p[i*m]);

    return 0;
}

int free2darr(float ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

/* TODO delete kai authn */
void DUMMYDUMDUM(int nx, int ny, float *u) {
int ix, iy;
int n=0;

for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*ny+iy) = n++;
}
