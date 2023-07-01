#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<mpi.h>

using namespace std;

int main() {
    int size, rank, i, j;
    int neighBor[1];
    neighBor[0] = MPI_PROC_NULL;
    //neighBor[1] = MPI_PROC_NULL;
    MPI_Comm comm;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm ,&size);
    MPI_Comm_rank(comm ,&rank);
    MPI_Datatype column;
    MPI_Status status;
    MPI_Request request;
    double test[2][6] = {{0.0, 1.0, 2.0, 4.0, 3000.0, 0.0},{0.0, 5.0, 6.0, 9.0, 7000.0, 0.0}};
    //double test1[2][4] = {{0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0}};

    int *xs = new int[size];
    int *xe = new int[size];
    int *ys = new int[size];
    int *ye = new int[size];
    
    xs[0] = 0;
    xe[0] = 1;

    xs[1] = 0;
    xe[1] = 1;

    ys[0] = 1;
    ye[0] = 1;

    ys[1] = 4;
    ye[1] = 4;
    
    /*if (rank == 1){
        for (i = xs[rank]; i <= xe[rank]; i++){
            for (j = ys[rank]; j <= ye[rank]; j++)
                cout << test[i][j] << " ";
            cout << endl;
        }
    }*/
    for (i = 0; i < 2; i++){
            for (j = 0; j < 6; j++)
                cout << test[i][j] << " ";
            cout << endl;
    }
    MPI_Type_vector(2, 1, 6, MPI_DOUBLE , &column);
    MPI_Type_commit(&column);
    int flag = 1, flag2 = 2;
    /*if (rank == 0)
        MPI_Sendrecv(&test[xs[rank]][ye[rank]], 1, column, 1, flag , &test[xs[rank]][ys[rank]-1], 1, column, 1, flag, comm, &status);
    if (rank == 1)
        MPI_Sendrecv(&test[xs[rank]][ys[rank]], 1, column, 0, flag , &test[xs[rank]][ye[rank]+1], 1, column, 0, flag, comm, &status);*/
    if (rank == 1)
        MPI_Send(&test[xs[1]][ys[1]], 1, column, 0, flag, comm);
    if (rank == 0)
        MPI_Recv(&test[xs[0]][ye[0]+1], 1, column, 1 , flag,  comm , &status);
    //cout << endl;
    if (rank == 0){
        cout << /*test[xs[rank]][ye[rank]] =*/ test[xe[rank]][ye[rank]+1] << test[xs[rank]][ye[rank]+1] << endl;
        for (i = 0; i < 2; i++){
                for (j = 0; j < 6; j++)
                    cout << test[i][j] << " ";
                cout << endl;
        }
    }
    delete[] xs;
    delete[] ys;
    delete[] xe;
    delete[] ye;
    MPI_Type_free(&column);
    MPI_Finalize();
    return 0;
}