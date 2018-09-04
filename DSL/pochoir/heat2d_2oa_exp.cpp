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

/* Test bench - 2D heat equation, Non-periodic version */
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

#define CHECK_RESULTS 1

/**
 *  * Get current time in seconds.
 *   */
double seconds ()
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  return ((double) tv.tv_sec) + 1e-6 * tv.tv_usec;
}

void check_result(int t, int i, int j, double a, double b)
{
	if (abs(a - b) < TOLERANCE) {
//		printf("a(%d, %d) == b(%d, %d) == %f : passed!\n", t, i, t, i, a);
	} else {
		printf("a(%d, %d, %d) = %f, b(%d, %d, %d) = %f : FAILED!\n", t, i,j, a, t, i,j, b);
	}
}

#define mod(r, m) (( r) %( m) + (( r) <0)? (m) :0)

Pochoir_Boundary_2D(heat_bv_2D, a, t, i,j)
  double val;

  if(t%2 == 0) {
    val = 1. + i*0.1 + j*0.01;
  } else {
    val = 2. + i*0.1 + j*0.01;
  }

  return val;

    //return 0;
Pochoir_Boundary_End

void print_at_timestep(int t, int i, int j, Pochoir_Array_2D(double) *arr) {
  printf("Arr[%d;%d,%d] = %f\n", t, i, j, arr->interior(t,i,j));
}

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
    N_SIZE = StrToInt(argv[1])-1;
    M_SIZE = StrToInt(argv[2])-1;
    O_SIZE = StrToInt(argv[3]);
    T_SIZE = StrToInt(argv[4]);
    printf("N_SIZE = %d, M_SIZE = %d, O_SIZE = %d, T_SIZE = %d\n", N_SIZE, M_SIZE, O_SIZE, T_SIZE);
	/* data structure of Pochoir - row major */
    Pochoir_Shape_2D heat_shape_2D[] = {{1, 0, 0}, {0, 1, 0}, {0,-1,0}, {0,0,1}, {0, 0, -1}, {0,0,0}};

    Pochoir_Array_2D(double) a(N_SIZE,M_SIZE);

#if CHECK_RESULTS
    Pochoir_Array_2D(double) b(N_SIZE,M_SIZE);
    //Pochoir_Array_2D(double) b(N_SIZE+2,M_SIZE+2);
    b.Register_Shape(heat_shape_2D);
#endif

    Pochoir_2D heat_2D(heat_shape_2D);

    Pochoir_Kernel_2D(heat_2D_fn, t, i,j)
	   a(t+1, i, j) = 0.125 * (
         a(t, i+1, j) -2.0*a(t,i,j) + a(t, i-1,j) 
         + a(t,i,j-1) - 2.0*a(t,i,j) + a(t,i,j+1)) + a(t,i,j);
    Pochoir_Kernel_End

    a.Register_Boundary(heat_bv_2D);
    heat_2D.Register_Array(a);

	for (int i = 0; i < N_SIZE; ++i) {
    for(int j = 0; j < M_SIZE; ++j) {
      a(0,i,j) = 1. + (i*0.1) + (j*0.01);
      a(1,i,j) = 2. + (i*0.1) + (j*0.01);
#if CHECK_RESULTS
      //b(0,i+1,j+1) = a(0,i,j);
      //b(1,i+1,j+1) = a(1,i,j);
      b(0,i,j) = a(0,i,j);
      b(1,i,j) = a(1,i,j);
#endif
    }
	} 
  double t1, t2;

#if 1
  t1 = seconds();
  heat_2D.Run(T_SIZE, heat_2D_fn);
  t2 = seconds();

    double nflops = (double) (N_SIZE - 2) * (double) (M_SIZE - 2) * T_SIZE * 9.0;
  cout << "FLOPs in stencil code: " << nflops << endl;
	cout << "Time spent in stencil coe: " << t2-t1 << " s" << endl;
  cout << "Performance in GFLOP/s: " << nflops / (1e9 * (t2-t1)) << endl;

#endif

  /*
   * Check results
   */
#if CHECK_RESULTS
  printf("Doing Naive Loop\n");
  for (int times = 0; times < TIMES; ++times) {

    for (int t = 0; t < T_SIZE; ++t) {
      for(int i = 1; i < N_SIZE-1; ++i) {
        for (int j = 1; j < M_SIZE-1; ++j) {

          b.interior(t+1, i, j) = 0.125 * (
              b.interior(t, i+1, j) -2.0*b.interior(t,i,j) + b.interior(t, i-1,j) 
              + b.interior(t,i,j-1) - 2.0*b.interior(t,i,j) + b.interior(t,i,j+1)) + b.interior(t,i,j);

        } 
      }
    }
  }

//  printf("Checking Results\n");
  t = T_SIZE - 1;
//  for (int i = 0; i < N_SIZE; ++i) {
//    for (int j = 0; j < M_SIZE; ++j) {
//      check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
//    } 
//  }

  //printf("a.interior(%d,%d,%d)=%f\n",t,2,2,a.interior(t,2,2));
  print_at_timestep(t,1,1,&a);
  print_at_timestep(t,2,2,&a);
  print_at_timestep(t,3,3,&a);
  //print_at_timestep(t,2+1,2+1,&b);
  print_at_timestep(t,3,3,&b);
#endif

	return 0;
}
