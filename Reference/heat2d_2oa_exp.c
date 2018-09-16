/*
 * Reference for heat two-dimensional 2nd order accurate (spatial) explicit method
 *
 * Original equation: U_tt = U_xx + U_yy
 * Solved by: u(t+1,x,y) = u(t,x,y)
 *  + 0.125 * (
 *    u(t,x-1,y) - 2(t,x,y) + u(t,x+1,y)
 *    u(t,x,y-1) - 2(t,x,y) + u(t,x,y+1)
 *  )
 *
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 8-26-2018
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define IDX(i,j) ((i)*x_max + (j))

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

    alpha = 1.f / (float) x_max;
    beta = 2.f / (float) y_max;

    /* initialize the first timesteps */
	#pragma omp parallel for private (k,j,i)
    for (j = 0; j < y_max; j++)
    {
      for (i = 0; i < x_max; i++)
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
      for (y = 1; y < y_max - 1; y++)
      {
        for (x = 1; x < x_max - 1; x++)
        {
          u_0_1[IDX(x, y)] = u_0_0[IDX(x, y)] +
            0.125 * (
                u_0_0[IDX(x+1, y)] -2*u_0_0[IDX(x,y)] + u_0_0[IDX(x-1, y)] +
                u_0_0[IDX(x, y+1)] -2*u_0_0[IDX(x,y)] + u_0_0[IDX(x, y-1)]);

        }
      }

    	float* tmp = u_0_0;
    	u_0_0 = u_0_1;
    	u_0_1 = tmp;
	}
	t2 = seconds ();

    /* print statistics */    
    nFlops = (double) (x_max-2) * (double) (y_max-2) * T_MAX * 9.0;
    printf ("FLOPs in stencil code:      %e\n", nFlops);    
	printf ("Time spent in stencil code: %f\n", t2 - t1);
	printf ("Performance in GFlop/s:     %f\n", nFlops / (1e9 * (t2 -t1)));

  printf("Value is u_0_0[2,2]=%f, u_0_1[2,2]=%f\n", u_0_0[IDX(2,2)], u_0_1[IDX(2,2)]);
   
    /* clean up */
	free (u_0_0);
	free (u_0_1);
	
	return EXIT_SUCCESS;
}


