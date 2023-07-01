#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926535

using namespace std;

double F(double x) {
    double f;
    f = 4.0/(1 + pow(x,2));
    return f;
}

int main(int argc, char* argv[]){
    #ifdef _OPENMP
	    printf("OpenMP is supported! %d \n", _OPENMP);
    #endif     

    int N = 10, myid, num_threads;
    double sum = 0, dx, tmp, s=0;
    double x_start = 0.0, x_end = 1.0;
    double time_init, time_final, elapsed_time;
    double *x = new double[N];
    
    dx = (x_end - x_start)/(N-1);

    /*#pragma omp parallel shared(dx, x_start) private(myid, s)
    {   
        time_init = omp_get_wtime();
        myid = omp_get_thread_num();
        num_threads = omp_get_num_threads();

        #pragma omp for
            for (int i = 0; i < N-1; i++){
                s +=  (F(i * dx + x_start) + F((i+1) * dx + x_start))/2.0 * dx;
            }
        #pragma omp critical
            {
                sum = sum + s;
            }
        #pragma omp barrier
        time_final = omp_get_wtime();
    }*/

    #pragma omp parallel reduction(+:sum)
    {
        time_init = omp_get_wtime();
        myid = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        #pragma omp for
            for (int i = 0; i < N-1; i++){
                sum +=  (F(i * dx + x_start) + F((i+1) * dx + x_start))/2.0 * dx;
            }
        time_final = omp_get_wtime();
    }

    printf("Resulting integral = %f, number of threads = %d\n", sum, num_threads);
    printf("Elapsed time = %f\n", time_final - time_init);


    free(x);
    return 0;
}
