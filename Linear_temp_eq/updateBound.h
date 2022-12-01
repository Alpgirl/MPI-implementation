#include "mpi.h"

/*******************************************************************/
/*             Update Bounds of subdomain with rank process          */
/*******************************************************************/

void updateBound(double** x, int neighBor[], MPI_Comm comm2d, int rank, 
    int* xs, int* ys, int* xe, int* ye, int ycell) {
    int S = 0, E = 1, N = 2, W = 3;
    int flag;
    MPI_Status status; // program read field values MPI_SOURCE and MPI_TAG to identify process

    /****************** North/South communication ******************/
    flag = 1;
    /* Send my boundary to North and receive from South */
    MPI_Sendrecv(&x[xe[rank]][ys[rank]], ycell, MPI_DOUBLE, neighBor[N], flag, &x[xs[rank]-1][ys[rank]], ycell, 
                    MPI_DOUBLE, neighBor[S], flag, comm2d, &status);

    /* Send my boundary to South and receive from North */
    MPI_Sendrecv(&x[xs[rank]][ys[rank]], ycell, MPI_DOUBLE, neighBor[S], flag, &x[xe[rank]+1][ys[rank]], ycell, 
                    MPI_DOUBLE, neighBor[N], flag, comm2d, &status);

    /****************** East/West communication ********************/
    flag = 2;
    /* Send my boundary to East and receive from West */
    for (int i = xs[rank]-1; i <= xe[rank]+1; i++){
        MPI_Sendrecv(&x[i][ye[rank]], 1, MPI_DOUBLE, neighBor[E], flag, &x[i][ys[rank]-1], 1, MPI_DOUBLE, 
                    neighBor[W], flag, comm2d, &status);

    /* Send my boundary to West and receive from East */
        MPI_Sendrecv(&x[i][ys[rank]], 1, MPI_DOUBLE, neighBor[W], flag, &x[i][ye[rank]+1], 1, MPI_DOUBLE, 
                    neighBor[E], flag, comm2d, &status);
    }
}

void computeNext(double** x0, double** x, double dt, double h_x, double h_y, 
    int rank, int *xs, int* ys, int* xe, int* ye, float a){

    int i, j;
    double a1, a2, a3, a4;
    
    for (i = xs[rank]; i <= xe[rank]; i++){
        for(j = ys[rank]; j <= ye[rank]; j++){
            a1 = x0[i+1][j] - x0[i][j];
            a2 = - x0[i][j] + x0[i-1][j];
            a3 = x0[i][j+1] - x0[i][j];
            a4 = - x0[i][j] + x0[i][j-1];
            //if (i == xs[rank]+2 && j == ys[rank]+2) std::cout << x0[i+1][j] << " "  << x0[i][j] << " "<< x0[i-1][j] << " " << a1+a2<< std::endl;
            x[i][j] = dt * pow(a,2) * ((a1+a2)/pow(h_x,2) + (a3+a4)/pow(h_y,2)) + x0[i][j];
            //if (i == xs[rank]+2 && j == ys[rank]+2) std::cout << x0[i+1][j] << " "  << x0[i][j] << std::endl;
        }
    }
    for (i = xs[rank]; i <= xe[rank]; i++){
        for(j = ys[rank]; j <= ye[rank]; j++){
            x0[i][j] = x[i][j];
            
        }
    }
}