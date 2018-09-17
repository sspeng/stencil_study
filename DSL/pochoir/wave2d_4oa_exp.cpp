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

/* Test bench - 2D wave equation, 4th order accurate; explicit */
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

Pochoir_Boundary_2D(bv, arr, t, i,j)
  double val;

  if(t%2 == 0) {
    val = 1. + i*0.1 + j*0.01;
  } else {
    val = 2. + i*0.1 + j*0.01;
  }

  return val;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    double min_tdiff = INF;
    int N_SIZE = 0, M_SIZE, O_SIZE, T_SIZE = 0;

    if (argc < 5) {
        printf("argc < 5, quit! <x_max> <y_max> <z_max> <t_max> \n");
        exit(1);
    }
    N_SIZE = StrToInt(argv[1]);
    M_SIZE = StrToInt(argv[2]);
    O_SIZE = StrToInt(argv[3]);
    T_SIZE = StrToInt(argv[4]);
    printf("N_SIZE = %d, M_SIZE = %d, O_SIZE = %d, T_SIZE = %d\n", N_SIZE, M_SIZE, O_SIZE, T_SIZE);

	/* data structure of Pochoir - row major */
    Pochoir_Shape_2D shape_2D[] = {{1, 0, 0}, {0, 1, 0}, {0,-1,0}, {0,0,1}, {0, 0, -1}, {0,0,0}, {-1,0,0}};
    Pochoir_Array_2D(double) a(N_SIZE * M_SIZE * T_SIZE);

    Pochoir_2D wave_2D(shape_2D);

    double sc1 = 1.0/12.0;
    double sc2 = 4.0/3.0;
    double sc3 = 5.0/2.0;

    Pochoir_Kernel_2D(fn, t, i, j)
      a(t+1,i,j) = (
          (-sc1*a(t,i+2,j) + sc2*a(t,i+1,j) - sc3*a(t,i,j) + sc2*a(t,i-1,j) - sc1*a(t,i-2,j)) +
          (-sc1*a(t,i,j+2) + sc2*a(t,i,j+1) - sc3*a(t,i,j) + sc2*a(t,i,j-1) - sc1*a(t,i,j-2)))
        + (2*a(t,i,j)) - a(t,i,j);
    Pochoir_Kernel_End

    a.Register_Boundary(bv);
    wave_2D.Register_Array(a);

    for (int i = 0; i < N_SIZE; ++i) {
      for (int j = 0; j < M_SIZE; ++j) {
        a(0,i) = 1. + i*0.1 + j*0.01;
        a(1,i) = 2. + i*0.1 + j*0.01;
      }
	} 
  double t1, t2;

  t1 = seconds();
  wave_2D.Run(T_SIZE, fn);
  t2 = seconds();

    double nflops = (double) (N_SIZE - 2) * (double) (M_SIZE - 2) * T_SIZE * 22.0;
  cout << "FLOPs in stencil code: " << nflops << endl;
	cout << "Time spent in stencil code: " << t2-t1 << " s" << endl;
  cout << "Performance in GFLOP/s: " << nflops / (1e9 * (t2-t1)) << endl;

	return 0;
}
