/*
 * Reference for laplacian one-dimensional 4th order accurate (spatial) explicit method
 *
 * Original equation: U_xx = 0
 * Solved by: u(t+1,x) = 
 *    (1/12u(t,x-2) - 4/3u(t,x-1) + 5/2(t,x) - 4/3u(t,x+1) + 1/12u(t,x+2))
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 8-26-2018
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//#define IDX(i,j,k) ((i)+x_max*((j)+y_max*(k)))
#define IDX(i) (i)

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
    u_0_0 = (float*) malloc (x_max * sizeof (float));
    u_0_1 = (float*) malloc (x_max * sizeof (float));

    alpha = 1.f / (float) x_max;
    beta = 2.f / (float) x_max;
    //beta = 2.f / (float) y_max;

    /* initialize the first timesteps */
	#pragma omp parallel for private (k,j,i)
    for (i = 0; i < x_max; i++)
    {
      u_0_0[IDX(i)] = 1. + i*0.1;
      u_0_1[IDX(i)] = 2. + i*0.1;
    }

    double sc1 = 1.0/12;
    double sc2 = 4.0/3.0;
    double sc3 = 5.0/2.0;

    /* do the calculation */ 
	t1 = seconds();
	for (t = 0; t < T_MAX; t++)
	{
		#pragma omp parallel for private(z,y,x)
    for (x = 2; x < x_max - 2; x++)
    {

      u_0_1[IDX(x)] = (sc1*u_0_0[IDX(x-2)] -  sc2*u_0_0[IDX(x-1)] + sc3* u_0_0[IDX(x)] - sc2*u_0_0[IDX(x+1)] + sc1*u_0_0[IDX(x+2)]);
    }

    	float* tmp = u_0_0;
    	u_0_0 = u_0_1;
    	u_0_1 = tmp;
	}
	t2 = seconds ();

    /* print statistics */    
    nFlops = (double) (x_max-4) * T_MAX * 9.0;
    printf ("FLOPs in stencil code:      %e\n", nFlops);    
	printf ("Time spent in stencil code: %f\n", t2 - t1);
	printf ("Performance in GFlop/s:     %f\n", nFlops / (1e9 * (t2 -t1)));
   
    /* clean up */
	free (u_0_0);
	free (u_0_1);
	
	return EXIT_SUCCESS;
}


