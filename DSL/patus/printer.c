
    void print_at_timestep(int t, int i, int j, double *arr, int x_max) {
      int idx = i*x_max + j;
      printf("Arr[%d;%d,%d] (%d) = %f\n", t, i, j, idx, arr[idx]);
    }
    print_at_timestep(99, 2, 2, A_0_1, x_max);

