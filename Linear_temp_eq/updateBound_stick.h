#include "mpi.h"

/*******************************************************************/
/*             Update Bounds of subdomain with rank process          */
/*******************************************************************/

void updateBound2(double** x, int N_x_global,  MPI_Comm comm, int rank,int size,
    int* xs, int* ys, int* xe, int* ye) {
    int S = 0, E = 1, N = 2, W = 3;
    int flag = 0;
    MPI_Status status; // program read field values MPI_SOURCE and MPI_TAG to identify process
    if(size > 1){
        for (int i = xs[rank]; i <= xe[rank]; i++){
            if(rank==0) 
                    MPI_Sendrecv(&x[i][ye[rank]], 1, MPI_DOUBLE, rank+1, flag, &x[i][ye[rank]+1], 1, MPI_DOUBLE, rank+1, flag, comm , &status);
            else if (rank == size - 1) 
                    MPI_Sendrecv(&x[i][ys[rank]], 1, MPI_DOUBLE, rank-1, flag, &x[i][ys[rank]-1], 1, MPI_DOUBLE, rank-1, flag , comm , &status);
            else{
                MPI_Sendrecv(&x[i][ys[rank]], 1, MPI_DOUBLE, rank-1, flag, &x[i][ys[rank]-1], 1, MPI_DOUBLE, rank-1, flag , comm , &status);
                MPI_Sendrecv(&x[i][ye[rank]], 1, MPI_DOUBLE, rank+1, flag, &x[i][ye[rank]+1], 1, MPI_DOUBLE, rank+1, flag, comm , &status);
            }
        }
    }
}

void computeNext(double** x0, double** x, double dt, double h_x, double h_y, 
    int rank, int *xs, int* ys, int* xe, int* ye, double a){

    int i, j;
    double a1, a2, a3, a4;
    
    for (i = xs[rank]; i <= xe[rank]; i++){
        for(j = ys[rank]; j <= ye[rank]; j++){
            a1 = x0[i+1][j] - x0[i][j];
            a2 = - x0[i][j] + x0[i-1][j];
            a3 = x0[i][j+1] - x0[i][j];
            a4 = - x0[i][j] + x0[i][j-1];
            x[i][j] = dt * a * ((a1+a2)/pow(h_x,2) + (a3+a4)/pow(h_y,2)) + x0[i][j];
        }
    }
    for (i = xs[rank]; i <= xe[rank]; i++){
        for(j = ys[rank]; j <= ye[rank]; j++){
            x0[i][j] = x[i][j];
            
        }
    }
}