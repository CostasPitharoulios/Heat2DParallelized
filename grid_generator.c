#include <stdio.h>
#include <stdlib.h>

#define NXPROB      80                 /* x dimension of problem grid */
#define NYPROB      64                 /* y dimension of problem grid */

int main(int argc, char *argv[]){
    void inidat(), prtdat();
    float  u[2][NXPROB][NYPROB];        /* array for grid */
    int	taskid,                     /* this task's unique id */
	numworkers,                 /* number of worker processes */
	numtasks,                   /* number of tasks */
	rc;

    /* Initialize grid */
    printf("Initializing grid and writing initial.dat file...\n");
    inidat(NXPROB, NYPROB, u);
    prtdat(NXPROB, NYPROB, u, "initial.dat");
}

void inidat(int nx, int ny, float *u) {
    int ix, iy;

    for (ix = 0; ix <= nx-1; ix++) 
      for (iy = 0; iy <= ny-1; iy++)
         *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

void prtdat(int nx, int ny, float *u1, char *fnam) {
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "wb");

    fwrite(u1, sizeof(float), ny*nx, fp);

    fclose(fp);
}
