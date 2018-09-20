/*
 * Reference for laplacian two-dimensional 2nd order accurate (spatial) explicit method
 *
 * Original equation: U_xx + U_yy = 0
 * Solved by: u(t+1,x,y) = 
 *      (u(t,x,y) - u(t,x-1,y) - u(t,x+1,y) - u(t,x,y-1) u(t,x,y+1))/4.0
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 8-26-2018
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define IDX(i,j) ((int) (((int) i)*x_max) + ((int) j))

#if defined(_OPENMP)
#	include <omp.h>
#endif


/**
 * Get current time in seconds.
 */
double seconds ()
{
    struct timeval tv;
    gettimeofday (&tv, NULL);
    return ((double) tv.tv_sec) + 1e-6 * tv.tv_usec;
}

/**
 * Do the calculation.
 */
int main(int argc, char** argv)
{
    int x_max, y_max, z_max;
    int i, j, k, t;
    int x, y, z;
    int T_MAX;
    double t1, t2, nFlops;

    float* __restrict__ u_0_0 = NULL;
    float* __restrict__ u_0_1 = NULL;
    float alpha, beta;

	if (argc != 5)
	{
		printf ("Wrong number of parameters.\n", argv[0]);
		exit (-1);
	}
	
	x_max = atoi (argv[1]);
	y_max = atoi (argv[2]);
	z_max = atoi (argv[3]);
	T_MAX = atoi (argv[4]);

 
    /* allocate memory */
    u_0_0 = (float*) malloc (x_max * y_max * sizeof (float));
    u_0_1 = (float*) malloc (x_max * y_max * sizeof (float));


    /* initialize the first timesteps */
	#pragma omp parallel for private (k,j,i)
    for (i = 0; i < x_max; i++)
    {
      for (j = 0; j < y_max; j++)
      {
        u_0_0[IDX(i,j)] = 1. + i*0.1 + j*0.01;
        u_0_1[IDX(i,j)] = 2. + i*0.1 + j*0.01;
      }
    }

    /* do the calculation */ 
	t1 = seconds();
	for (t = 0; t < T_MAX; t++)
	{
		#pragma omp parallel for private(z,y,x)
    for (x = 1; x < x_max - 1; x++)
    {
      for (y = 1; y < y_max - 1; y++)
      {
        u_0_1[IDX(x,y)] = (u_0_0[IDX(x,y)] 
            - u_0_0[IDX(x-1,y)] - u_0_0[IDX(x+1,y)]
            - u_0_0[IDX(x,y-1)] - u_0_0[IDX(x,y+1)])/4.0;
      }
    }

    	float* tmp = u_0_0;
    	u_0_0 = u_0_1;
    	u_0_1 = tmp;
	}
	t2 = seconds ();

    /* print statistics */    
    nFlops = (double) (x_max-2) * (double) (y_max-2) * T_MAX * 5.0;
    printf ("FLOPs in stencil code:      %e\n", nFlops);    
	printf ("Time spent in stencil code: %f\n", t2 - t1);
	printf ("Performance in GFlop/s:     %f\n", nFlops / (1e9 * (t2 -t1)));
   
    /* clean up */
	free (u_0_0);
	free (u_0_1);
	
	return EXIT_SUCCESS;
}


