#include "mpi.h"

/*******************************************************************/
/*             Update Bounds of subdomain with rank process          */
/*******************************************************************/

void updateBound(double** x, int neighBor[], MPI_Comm comm2d, MPI_Datatype column_type, int rank, 
    int* xs, int* ys, int* xe, int* ye, int ycell) {
    int S = 0, E = 1, N = 2, W = 3;
    int flag;
    MPI_Status status; // program read field values MPI_SOURCE and MPI_TAG to identify process

    /****************** North/South communication ******************/
    flag = 1;
    /* Send my boundary to North and receive from South */
    MPI_Sendrecv(&x[xe[rank]][ys[rank]], ycell, MPI_DOUBLE, neighBor[N], flag, &x[xs[rank]-1][ys[rank]], ycell, 
                    MPI_DOUBLE, neighBor[S], flag, comm2d, &status);
    //MPI_Sendrecv( const void* sendbuf , MPI_Count sendcount , MPI_Datatype sendtype , int dest , int sendtag , void* recvbuf , 
     //             MPI_Count recvcount , MPI_Datatype recvtype , int source , int recvtag , MPI_Comm comm , MPI_Status* status);

    /* Send my boundary to South and receive from North */
    MPI_Sendrecv(&x[xs[rank]][ys[rank]], ycell, MPI_DOUBLE, neighBor[S], flag, &x[xe[rank]+1][ys[rank]], ycell, 
                    MPI_DOUBLE, neighBor[N], flag, comm2d, &status);

    /****************** East/West communication ********************/
    flag = 2;
    /* Send my boundary to East and receive from West */
    for (int i = xs[rank]-1; i <= xe[rank]+1; i++){
        MPI_Sendrecv(&x[i][ye[rank]], 1, /*column_type*/ MPI_DOUBLE, neighBor[E], flag, &x[i][ys[rank]-1], 1, /*column_type*/MPI_DOUBLE, 
                    neighBor[W], flag, comm2d, &status);

    /* Send my boundary to West and receive from East */
        MPI_Sendrecv(&x[i][ys[rank]], 1, /*column_type*/MPI_DOUBLE, neighBor[W], flag, &x[i][ye[rank]+1], 1, /*column_type*/MPI_DOUBLE, 
                    neighBor[E], flag, comm2d, &status);
    }
    /*if (rank == 1)
        MPI_Send(&x[xs[1]][ys[1]], 1, column_type, 0, flag, comm2d);*/
        /*for (int i = xs[1]; i <= xe[1]; i++)
            std::cout << x[i][ys[1]-1] << " ";
        std::cout << std :: endl;*/
    /*if (rank == 0)
        MPI_Recv(&x[xs[0]][ye[0]+1], 1, column_type, 1 , flag,  comm2d , &status);*/
}

void computeNext(double** x0, double** x, double dt, double h_x, double h_y, double* diff, 
    int rank, int *xs, int* ys, int* xe, int* ye, double a){

    int i, j;
    double ldiff;
    
    for (i = xs[rank]; i <= xe[rank]; i++){
        for(j = ys[rank]; j <= ye[rank]; j++){
            //if (i == xs[rank] and j == ys[rank]+1) std::cout << x0[i][j] << " " << x0[i-1][j] << " " << x0[i][j-1] << " " << x0[i+1][j] << " " << x0[i][j+1] << std ::endl;
            x[i][j] = pow(a,2)*dt/pow(h_x,2) * (x0[i+1][j] - 2.*x0[i][j] + x0[i-1][j]) + pow(a,2)*dt/pow(h_y,2) * (x0[i][j+1] - 2.*x0[i][j] + x0[i][j-1]) 
                + x0[i][j];
           // if (i == xs[rank] and j == ys[rank]+1) std::cout << x[i][j] << std :: endl;
        }
    }
    *diff = 0.0;
    for (i = xs[rank]; i <= xe[rank]; i++){
        for(j = ys[rank]; j <= ye[rank]; j++){
            ldiff = x0[i][j] - x[i][j];
            *diff += ldiff * ldiff;
            x0[i][j] = x[i][j];
            //if (i == xs[rank] and j == ys[rank]+1) std::cout << x0[i][j] << std :: endl;
            
        }
    }
}