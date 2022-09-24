//
//  hello_w_ABC.cpp
//  
//
//  Created by Inna on 24/09/2022.
//

#include "hello_w_ABC.hpp"
#include <mpi.h>
#include <iostream>
#include <string.h>
using namespace std;

int main() {
    MPI_Init(NULL, NULL); /* starts MPI */
    int rank;
    int size;
    const char * names[3] = {"Alice", "Bob", "Kate"};
    const char * greeting = "Hello, ";
    char message[15];
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get current process id */
    // MPI_COMM_WORLD is constructed my MPI. It encloses all of the processes in the job and returns the amount of processe thet requested for the job
    
    MPI_Comm_size(MPI_COMM_WORLD, &size); /* get number of processes */
    
    if (rank != 0) {
        strcpy(message, greeting);
        strcat(message, names[rank]);
        // Note that after MPI_Recv address of message changes!
        cout << names[rank] << " received message '" << message << "' from " << names[rank-1] << endl;
        
        MPI_Recv(&message, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    MPI_Send(&message, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        strcpy(message, greeting);
        strcat(message, names[rank]);
        
        cout << names[rank] << " received message '" << message << "' from " << names[size-1] << endl;
        
        MPI_Recv(&message, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    MPI_Finalize();
    return 0;
}

