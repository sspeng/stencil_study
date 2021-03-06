/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 * 		                     Charles E. Leiserson <cel@mit.edu>
 * 	 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* Test bench - 3D heat equation, Non-periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define N_RANK 2
#define TOLERANCE (1e-6)

/**
 *  * Get current time in seconds.
 *   */
double seconds ()
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  return ((double) tv.tv_sec) + 1e-6 * tv.tv_usec;
}

void check_result(int t, int i, double a, double b)
{
	if (abs(a - b) < TOLERANCE) {
//		printf("a(%d, %d) == b(%d, %d) == %f : passed!\n", t, i, t, i, a);
	} else {
		printf("a(%d, %d) = %f, b(%d, %d) = %f : FAILED!\n", t, i, a, t, i, b);
	}

}

Pochoir_Boundary_3D(heat_bv_3D, arr, t, i,j,k)
  double val;

  if(t%2 == 0) {
    val = 1. + i*0.1 + j*0.01 + k*0.001;
  } else {
    val = 2. + i*0.1 + j*0.01 + k*0.001;
  }

  return val;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    double min_tdiff = INF;
    int N_SIZE = 0, T_SIZE = 0, M_SIZE=0, O_SIZE=0;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    N_SIZE = StrToInt(argv[1]) - 1;
    M_SIZE = StrToInt(argv[2]) - 1;
    O_SIZE = StrToInt(argv[3]) - 1;
    T_SIZE = StrToInt(argv[4]);
    printf("N_SIZE = %d, M_SIZE = %d, O_SIZE = %d, T_SIZE = %d\n", N_SIZE, M_SIZE, O_SIZE, T_SIZE);
	/* data structure of Pochoir - row major */
    Pochoir_Shape_3D heat_shape_3D[] = {{1, 0, 0, 0}, {0, 1, 0,0}, {0,-1,0,0}, {0,0,1,0}, {0, 0, -1,0}, {0,0,0,1},{0,0,0,-1},{0,0,0,0}};
    Pochoir_Array_3D(double) a(N_SIZE,M_SIZE,O_SIZE);
    Pochoir_3D heat_3D(heat_shape_3D);

    Pochoir_Kernel_3D(heat_3D_fn, t, i,j,k)
	   a(t+1, i, j,k) = 0.125 * (
         a(t, i+1, j,k) -2.0*a(t,i,j,k) + a(t, i-1,j,k) 
         + a(t,i,j-1,k) -2.0*a(t,i,j,k) + a(t,i,j+1,k)
         + a(t,i,j,k-1) -2.0*a(t,i,j,k) + a(t,i,j,k+1)
         )
     + a(t,i,j,k);
    Pochoir_Kernel_End

    a.Register_Boundary(heat_bv_3D);
    heat_3D.Register_Array(a);

    for (int i = 0; i < N_SIZE; ++i) {
      for(int j = 0; j < M_SIZE; ++j) {
        for(int k = 0; k < O_SIZE; ++k) {
          a(0,i,j,k) = 1. + i*0.1 + j*0.01 + k*0.001;
          a(1,i,j,k) = 2. + i*0.1 + j*0.01 + k*0.001;
        }
      }
    } 
    double t1, t2;

#if 1
  t1 = seconds();
  heat_3D.Run(T_SIZE, heat_3D_fn);
  t2 = seconds();

    double nflops = (double) (N_SIZE - 2) * (double) (M_SIZE - 2) * (double) (O_SIZE - 2) * T_SIZE * 13.0;
  cout << "FLOPs in stencil code: " << nflops << endl;
	cout << "Time spent in stencil coe: " << t2-t1 << " s" << endl;
  cout << "Performance in GFLOP/s: " << nflops / (1e9 * (t2-t1)) << endl;

#endif

	return 0;
}
