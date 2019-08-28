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
 *   process and then distributed by blocks to the worker processes. At each 
 *   time step, worker processes must exchange border data with neighbors, 
 *   because a grid point's current temperature depends upon it's previous
 *   time step value plus the values of the neighboring grid points. Upon
 *   completion of all time steps, the worker processes return their results
 *   to the master process.
 *
 *   Two data files are produced: an initial data set and a final data set.
 * AUTHORS: Costas Pitharoulios, Simon Iyamu
 *   
 ****************************************************************************/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NXPROB      8                  /* x dimension of problem grid */
#define NYPROB      12                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define BEGIN       2                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};

void inidat(), prtdat(), updateExternal(), updateInternal(),  myprint(), DUMMYDUMDUM();
int malloc2darr(),free2darr(),isPrime(),checkSize();

int main (int argc, char *argv[]){
    float **local[2];               /* stores the block assigned to current task, surrounded by halo points */
    int	taskid,                     /* this task's unique id */
        numworkers,                 /* number of worker processes */
        dest, source,               /* to - from for message send-receive */
        left,right,up,down,         /* neighbor tasks */
        msgtype,                    /* for message types */
        xdim, ydim,                 /* dimensions of grid partition (e.x. 4x4) */
        rows, columns,              /* number of rows/columns of each block (e.x. 20x12) */
        i,j,x,y,ix,iy,iz,it;        /* loop variables */
    double start,finish;
    char inputfile[80] = "initial.dat";
    char outputfile[80] = "final.dat";
    MPI_Status status;

    /* Read arguments */
    for(i=1; i<argc; i++){
        if(!strcmp(argv[i],"-i"))
            strcpy(inputfile,argv[i+1]);
        if(!strcmp(argv[i],"-o"))
            strcpy(outputfile,argv[i+1]);
    }

    /* First, find out my taskid and how many tasks are running */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numworkers);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers;

    if (taskid == MASTER) {
        /************************* Master code *******************************/

        if ((isPrime(numworkers))){
            printf("ERROR: the number of workers is prime (%d).\n",numworkers);
            MPI_Abort(MPI_COMM_WORLD, 22);
            exit(22);
        }

        /* If the number of cells is not divisible by numworkers, abort */
        if ((NXPROB*NYPROB)%numworkers){
            printf("ERROR: number of cells is not divisible by the number of workers\n");
            MPI_Abort(MPI_COMM_WORLD, 22);
            exit(22);
        }

        if (!checkSize(inputfile)){
            printf("ERROR: grid size of input file is diffrent that the definitions of NXPROB, NYPROB or file doesn't exist\n");
            MPI_Abort(MPI_COMM_WORLD, 22);
            exit(22);
        }

        printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

        printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
#if 0
        for (ix=0; ix<NXPROB; ix++){
            for (j=0; j<NYPROB; j++)
                printf("%6.1f ", u[0][ix][j]);
            printf("\n\n");
        }
#endif

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
        for (i=1; i<numworkers; i++){

            /* Find the neighbours of this block */
            if (i < ydim) // if this is the first row
                up = MPI_PROC_NULL;
            else
                up = i - ydim;

            if (i >= ((xdim-1) * ydim)) //if this is the last row
               down = MPI_PROC_NULL;
            else
               down = i + ydim;

            if (i%ydim == 0) // if this is the first column
                left = MPI_PROC_NULL;
            else
                left = i-1;

            if (i%ydim == ydim-1)	//if this is the last column
                right = MPI_PROC_NULL;
            else
                right = i+1;

            /*  Now send startup information to each worker  */
            dest = i;
            MPI_Send(&xdim, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&ydim, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&columns, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&down, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        }

        /* Master computes its neighbours */
        left = MPI_PROC_NULL;
        up = MPI_PROC_NULL;
        if (numworkers == 1)
            right = down = MPI_PROC_NULL;
        else{
            right = 1;
            down = ydim;
        }

    }else{
        /* taskid != MASTER */

        /* Receive my offset, rows, neighbors and grid partition from master */
        source = MASTER; msgtype = BEGIN;
        MPI_Recv(&xdim, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&ydim, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&columns, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&down, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    }
    printf("LOG: Process %d: left:%d, right:%d, up:%d, down:%d\n",taskid,left,right,up,down);

    /* Define a new communicator with cartesian topology information, for communication optimization */
    MPI_Comm comm_cart;
    int dim[2] = {xdim,ydim}, period[2] = {0,0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 0, &comm_cart);

    /* Allocate contigious memory for the 2d arrays local[0] and local[1] */
    malloc2darr(&local[0], rows+2, columns+2);
    malloc2darr(&local[1], rows+2, columns+2);

    /* Initialize with 0's */
    for (iz=0; iz<2; iz++)
        for (ix=0; ix<rows+2; ix++) 
            for (iy=0; iy<columns+2; iy++) 
                local[iz][ix][iy] = 0.0;

    /* Preparing the datatypes for Parallel I/O */

    /* Define the datatype of send buffer elements */
    int sendsizes[2]    = {NXPROB, NYPROB};    /* grid size */
    int sendsubsizes[2] = {rows, columns};     /* local size without halo */
    int sendstarts[2]   = {0,0};

    MPI_Datatype type, sendsubarrtype;
    MPI_Type_create_subarray(2, sendsizes, sendsubsizes, sendstarts, MPI_ORDER_C, MPI_FLOAT, &type);
    MPI_Type_create_resized(type, 0, columns*sizeof(float), &sendsubarrtype);
    MPI_Type_commit(&sendsubarrtype);

    /* Define the datatype of receive buffer elements */
    int recvsizes[2]    = {rows+2, columns+2};         /* local array size */
    int recvsubsizes[2] = {rows, columns};          /* local size without halo */
    int recvstarts[2]   = {1,1};

    MPI_Datatype recvsubarrtype;
    MPI_Type_create_subarray(2, recvsizes, recvsubsizes, recvstarts, MPI_ORDER_C, MPI_FLOAT, &recvsubarrtype);
    MPI_Type_commit(&recvsubarrtype);
   
    /* Open file for reading */ 
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, inputfile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  
    /* Set view in order to define which portion of the file is visible to each worker */
    MPI_Offset disp = ((taskid/ydim)*rows*NYPROB + taskid%ydim*columns)*sizeof(float);
    MPI_File_set_view(fh, disp, MPI_FLOAT, sendsubarrtype, "native", MPI_INFO_NULL);

    /* Read from the file */
    MPI_File_read(fh, &(local[0][0][0]), 1, recvsubarrtype, &status);
    MPI_File_close(&fh);

    /// *** WORK STARTS HERE *** ///

    /* Start the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    iz = 0;

    MPI_Request RRequestR, RRequestL, RRequestU, RRequestD;
    MPI_Request SRequestR, SRequestL, SRequestU, SRequestD;

    /* Datatypes for matrix column */
    MPI_Datatype column; 
    MPI_Type_vector(rows, 1,columns+2, MPI_FLOAT, &column);
    MPI_Type_commit(&column);

    /* Requests for persistent communication */
    MPI_Request req[8];
    MPI_Status  stat[8];

    MPI_Recv_init(&(local[iz][1][0]), 1, column, left, 0, comm_cart, &(req[0]));
    MPI_Recv_init(&(local[iz][1][columns+1]), 1, column, right, 0, comm_cart, &(req[1]));
    MPI_Recv_init(&(local[iz][rows+1][1]), columns, MPI_FLOAT, down, 0, comm_cart, &(req[2])); 
    MPI_Recv_init(&(local[iz][0][1]), columns, MPI_FLOAT, up,0, comm_cart, &(req[3])); 

    MPI_Send_init(&(local[iz][1][columns]), 1, column, right, 0, comm_cart, &req[4]);
    MPI_Send_init(&(local[iz][1][1]), 1, column, left , 0, comm_cart, &req[5]);
    MPI_Send_init(&(local[iz][1][1]), columns, MPI_FLOAT, up, 0, comm_cart, &req[6]);
    MPI_Send_init(&(local[iz][rows][1]), columns, MPI_FLOAT, down ,0, comm_cart, &req[7]);
    
    MPI_Startall(8,req);
    MPI_Waitall(8,req,MPI_STATUS_IGNORE);

    for (it = 1; it <= STEPS; it++){

        /// *** RECEIVING PROCEDURES *** ///
        MPI_Irecv(&(local[iz][1][0]), 1, column, left, 0, comm_cart, &RRequestL); ///WARNING: 0??
        MPI_Irecv(&(local[iz][1][columns+1]), 1, column, right, 0, comm_cart, &RRequestR); ///WARNING: 0?
        MPI_Irecv(&(local[iz][rows+1][1]), columns, MPI_FLOAT, down, 0, comm_cart, &RRequestD); ///WARNING: 0??
        MPI_Irecv(&(local[iz][0][1]), columns, MPI_FLOAT, up,0, comm_cart, &RRequestU); ///WARNING: 0??

        /// *** SENDING PROCEDURES *** ///
        MPI_Isend(&(local[iz][1][columns]), 1, column, right, 0, comm_cart, &SRequestR);  //sends column to RIGHT neighbor
        MPI_Isend(&(local[iz][1][1]), 1, column, left , 0, comm_cart, &SRequestL);	//sends column to left neighbor
        MPI_Isend(&(local[iz][1][1]), columns, MPI_FLOAT, up, 0, comm_cart, &SRequestU);  //sends to UP neighbor
        MPI_Isend(&(local[iz][rows][1]), columns, MPI_FLOAT, down ,0, comm_cart, &SRequestD); //sends to DOWN neighbor

        /// *** CALCULATION OF INTERNAL DATA *** ///
        updateInternal(2, rows-1, columns,&local[iz][0][0], &local[1-iz][0][0]); // 2 and xdim-3 because we want to calculate only internal nodes of the block.
        //line 0 contains neighbor's values and line 1 is the extrnal line of the block, so we don't want them. The same for the one before last and the last line.

        if (right != MPI_PROC_NULL) MPI_Wait(&RRequestR , MPI_STATUS_IGNORE );
        if (left != MPI_PROC_NULL) MPI_Wait(&RRequestL , MPI_STATUS_IGNORE );
        if (up !=  MPI_PROC_NULL) MPI_Wait(&RRequestU , MPI_STATUS_IGNORE );
        if (down !=  MPI_PROC_NULL) MPI_Wait(&RRequestD , MPI_STATUS_IGNORE );

        /// *** CALCULATION OF EXTERNAL DATA *** ///
        updateExternal(1,rows, columns,right,left,up,down, &local[iz][0][0], &local[1-iz][0][0]);

        iz = 1-iz; 

        if (right != MPI_PROC_NULL) MPI_Wait(&SRequestR , MPI_STATUS_IGNORE );
        if (left != MPI_PROC_NULL) MPI_Wait(&SRequestL , MPI_STATUS_IGNORE );
        if (up !=  MPI_PROC_NULL) MPI_Wait(&SRequestU , MPI_STATUS_IGNORE );
        if (down !=  MPI_PROC_NULL) MPI_Wait(&SRequestD , MPI_STATUS_IGNORE );

#if 0
        for ( i=0; i<numworkers; i++){
            if (taskid == i){
                printf("=========== To kommati tou %d meta thn antallagh =========\n",i);
                for (ix=0; ix<rows+2; ix++){
                    for (j=0; j<columns+2; j++)
                        printf("%6.1f ", local[1-iz][ix][j]);
                    printf("\n\n");
                }
	        printf("=========== To kommati tou %d meta thn UPDATE =========\n",i);
                for (ix=0; ix<rows+2; ix++){
                    for (j=0; j<columns+2; j++)
                        printf("%6.1f ", local[iz][ix][j]);
                    printf("\n\n");
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif


    }

    /// *** WORK COMPLETE *** ///

    /* Stop the timer */
    finish = MPI_Wtime();

    /* Gather it all back */

    /* Each worker writes to its portion of the file */
    MPI_File_open(MPI_COMM_WORLD, outputfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  
    /* Set view in order to define which portion of the file is visible to each worker */
    MPI_File_set_view(fh, disp, MPI_FLOAT, sendsubarrtype, "native", MPI_INFO_NULL);

    /* Write to the file */
    MPI_File_write(fh, &(local[iz][0][0]), 1, recvsubarrtype, &status);
    MPI_File_close(&fh);

    printf("Process:%d, Elapsed time: %e secs\n",taskid,finish-start);

    /* Free malloc'd memory */
    free2darr(&local[0]);
    free2darr(&local[1]);

    MPI_Type_free(&type);
    MPI_Type_free(&sendsubarrtype);
    MPI_Type_free(&recvsubarrtype);
    MPI_Type_free(&column);

    for(i=0; i<8 ; i++)
        MPI_Request_free(&(req[i]));
    
    MPI_Finalize();
    return 0;
}

/**************************************************************************
 *  subroutine update
/// gets start = 2, end = xdim-1, ny = ydim = number of block columns without 
/// the two which keep LEFT AND RIGHT neighbors' values
 ****************************************************************************/
void updateInternal(int start, int end, int ny, float *u1, float *u2)
{

   int ix, iy;
   for (ix = start; ix <= end; ix++){ 
      for (iy = 2; iy <= ny-1; iy++){
         *(u2+ix*(ny+2)+iy) = *(u1+ix*(ny+2)+iy)  + 
                          parms.cx * (*(u1+(ix+1)*(ny+2)+iy) +
                          *(u1+(ix-1)*(ny+2)+iy) - 
                          2.0 * *(u1+ix*(ny+2)+iy)) +
                          parms.cy * (*(u1+ix*(ny+2)+iy+1) +
                         *(u1+ix*(ny+2)+iy-1) - 
                          2.0 * *(u1+ix*(ny+2)+iy));
       }
    }
}


/**************************************************************************
 *  subroutine updateExternal
///gets start = 1, end = xdim, ny= ydim = number of block columns without 
///the two which keep LEFT AND RIGHT neighbors' values
 ****************************************************************************/
void updateExternal(int start, int end, int ny,int right, int left,int up,int down, float *u1, float *u2)
{
    int ix, iy;
    ny+=2;
    end+=2;
	/// *** CALCULATING FIRST EXTERNAL ROW *** ///
    if (up != MPI_PROC_NULL) //this is because if the block haw not an up neighbor we shouldnt's caclulate halo
       ix = start;
    else
	    ix = start+1;

    if (left != MPI_PROC_NULL) //this is because if the block haw not a left neighbor we shouldnt's caclulate halo
        iy = 1;
    else
        iy = 2;

    int endny;
    if (right !=  MPI_PROC_NULL)
        endny = ny-3;
    else
        endny = ny-3;

    for (; iy <= endny; iy++) 
         *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                          parms.cx * (*(u1+(ix+1)*ny+iy) +
                          *(u1+(ix-1)*ny+iy) - 
                          2.0 * *(u1+ix*ny+iy)) +
                          parms.cy * (*(u1+ix*ny+iy+1) +
                         *(u1+ix*ny+iy-1) - 
                          2.0 * *(u1+ix*ny+iy));

	/// CALCULATING LAST EXTERNAL ROW *** ///
    if (down != MPI_PROC_NULL)
       ix =  end-2;
    else
	ix = end-3;
    if (left != MPI_PROC_NULL) //this is because if the block haw not a left neighbor we shouldnt's caclulate halo
        iy = 1;
    else
        iy = 2;
    if (right !=  MPI_PROC_NULL)
        endny = ny-2;
    else
        endny = ny-3;
    for (; iy <= endny; iy++) 
         *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                          parms.cx * (*(u1+(ix+1)*ny+iy) +
                          *(u1+(ix-1)*ny+iy) - 
                          2.0 * *(u1+ix*ny+iy)) +
                          parms.cy * (*(u1+ix*ny+iy+1) +
                         *(u1+ix*ny+iy-1) - 
                          2.0 * *(u1+ix*ny+iy));
	/// *** CALCULATING FIRST EXTERNAL COLUMN *** ///

    if (up != MPI_PROC_NULL) //this is because if the block haw not an up neighbor we shouldnt's caclulate halo
       ix = start;
    else
        ix = start+1;


    if (left != MPI_PROC_NULL) //this is because if the block haw not a left neighbor we shouldnt's caclulate halo
        iy = 1;
    else
        iy = 2;
    
    int endloop;
    if (down != MPI_PROC_NULL)
       endloop = end -2;
    else
       endloop = end -3; 

    for (; ix<endloop; ix++)
        *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                          parms.cx * (*(u1+(ix+1)*ny+iy) +
                          *(u1+(ix-1)*ny+iy) - 
                          2.0 * *(u1+ix*ny+iy)) +
                          parms.cy * (*(u1+ix*ny+iy+1) +
                         *(u1+ix*ny+iy-1) - 
                          2.0 * *(u1+ix*ny+iy)); 
 /// *** CALCULATING LAST EXTERNAL COLUMN *** ///

   if (up != MPI_PROC_NULL) //this is because if the block haw not an up neighbor we shouldnt's caclulate halo
       ix = start;
   else
       ix = start+1; 

   if (right != MPI_PROC_NULL)
        iy = ny -2;
    else 
        iy = ny-3;

    if (down != MPI_PROC_NULL)
       endloop = end -2; //the down right corner is calculated from row calculation, so we don't need to calculate again
    else 
       endloop = end -3; // the down right corner is calculated from row calculation, so we don't need to calculate again
    //printf("end =%d, endloop=%d\n\n", end, endloop);
    for (; ix<endloop; ix++)
       *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + parms.cx * (*(u1+(ix+1)*ny+iy) +
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
    if (n==1)
        return 0;
    if (n==2)
        return 1;
    if (n%2==0)
        return 0;
    for (i=3;i*i<=n;i+=2)
        if (n%i==0) 
            return 0;
    return 1;
}

/* Checks if grid size of given file is the same as NXPROB x NYPROB */
int checkSize(const char *filename){
    FILE *fp;
    fp = fopen(filename, "rb");
    if(!fp)
        return 0;
    fseek(fp, 0L, SEEK_END);
    long int sz = ftell(fp);
    fclose(fp);

    return sz == NYPROB*NXPROB*sizeof(float);
}
#if 0 
int checkSize(const char *filename){
    int x=0,y=0;
    FILE *fp;
    char *line = NULL, *token;
    size_t len = 0;
    ssize_t read;

    fp = fopen(filename, "r");
    if (fp == NULL)
        return 0;

    while ((read = getline(&line, &len, fp)) != -1) {
        y++;

        /* Count tokens */
        x = 0;
        token = strtok(line, " ");
        while (token != NULL) {
            x++;
            token = strtok(NULL, " ");
        }

        /* Check if x is right */
        if (x != NXPROB){
            fclose(fp);
            free(line);
            return 0;
        }
    }

    fclose(fp);
    if (line)
        free(line);
    
    /* Check y is right */
    if(y != NYPROB)
        return 0;
    return 1;
}
#endif 

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
