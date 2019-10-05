#include "mpi.h"
#include <omp.h>
#include <stdio.h>


int main(int argc, char *argv[]){
//#if 0
   int numworkers, taskid,provided;
    MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD,&numworkers);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);


    printf("taskid=%d\n", taskid);
//#endif
//   int i,j;
   #pragma omp parallel num_threads(4)
   {
       int i,j;
       for (i=0; i< 2; i++){
           #pragma omp single
	   {
	       printf("Inside single, thread:%d\n\n",  omp_get_thread_num());
	   }
	   #pragma omp for
	   for (j=0; j<3; j++){
	       printf("Double for [%d][%d] | THREAD: %d\n\n", i,j, omp_get_thread_num());
	   }
	   #pragma omp single
	   {
	       printf("Getting out of for with  i=%d, THREAD: %d\n\n", i, omp_get_thread_num());
	   }
       } 
   }

}





