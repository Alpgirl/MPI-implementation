#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926535

using namespace std;

double F(double x) {
    double f;
    f = 4.0/(1 + pow(x,2));
    return f;
}

void InsertArr(double *x, int ind, int N_total, int N_shift){
    int i;
    for (i = N_shift; i > ind; i--)
        x[i] = x[i-1];
}

void FieldX(double *x, double size, int N, double dx, int N_global){
    int cnt, i;
    for (i = 0; i <= N; i++){
        x[i] = i * dx;
    }
    cnt = 0;
    for (i = N/size; i < N; i += N/size){
        InsertArr(x, i+cnt, N_global, cnt + N_global);
        cnt += 1;
    }
}

int main(int argc, char* argv[]){
    int rank, size, N = 1000, N_global, N_add, N_eq, cnt = 0;
    double sum = 0, dx, tmp;
    double x_start = 0.0, x_end = 1.0;
    double int_total;
    double time_init, time_final, elapsed_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    N_add = N - N/size * size; // остаток от неравномерного деления отрезка
    N_eq = N/size * size;

    dx = (x_end - x_start)/N;
    N_global = N_eq + size;

    double *x = new double[N_global];
    double *x_local = new double[N_global/size];

    FieldX(x, size, N, dx, N_global);

    //MPI_Scatter(x, N_global/size, MPI_DOUBLE, x_local, N_global/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    time_init = MPI_Wtime();
    for (int i = rank * N_global/size; i < (rank + 1) * N_global/size; i++){
        x_local[cnt] = x[i];
        cnt++;
    }


    for (int i = 0; i < N_global/size - 1; i++){
        sum += ( F(x_local[i]) + F(x_local[i+1]))/2.0 * dx;
    }
    if (rank == 0){
        for (int i = N_eq; i < N; i++){
            sum += ( F(i * dx) + F((i+1)*dx))/2.0 * dx;
        }
    }

    MPI_Reduce(&sum, &int_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Ending time 
    time_final = MPI_Wtime();

    // Elapsed time
    elapsed_time = time_final - time_init;

    if (rank == 0) {
        cout << "Result of calculating = " << int_total << endl;
        cout << "Time = " << elapsed_time<< endl;
    }

    free(x);
    free(x_local);
    MPI_Finalize();
    return 0;
}