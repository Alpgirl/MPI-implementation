#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#define N_x  21
#define N_y  21
#define M  N_x * N_y
#define h  5.0
#define L_x  10.0
#define L_y  20.0
#define T_1  1000.0
#define T_2  300.0
using namespace std;

/*double calc_side (double *** T, int i, int j, double **alpha, int t, int h_x, int h_y, int tau) {
    float test, test1, test2, test3;
    test = alpha[0][1];
    test1 = alpha[1][1];
    test2 = alpha[1][0];
    test3 = alpha[0][0];
    T[i][j][t + 1] = tau * ((alpha[0][0] * (T[i+1][j][t] - T[i][j][t]) + alpha[0][1] * (- T[i][j][t] + T[i-1][j][t]))/pow(h_x,2) + 
    (alpha[1][0] * (T[i][j+1][t] - T[i][j][t]) + alpha[1][1] * (- T[i][j][t] + T[i][j-1][t]))/pow(h_y,2)) + T[i][j][t];
    T[i][j][t+1] = tau * ((alpha[0][0] + alpha[0][1])/pow(h_x,2) + (alpha[1][0] + alpha[1][1])/pow(h_y,2)) + T[i][j][t];
    return T[i][j][t+1];
}*/
void initMat(double *T, float h_x, float h_y){
    int i, j;
    for(i = 0; i < M; i++)
        for (j = 0; j < M; j++){
            if (h_x * i < h && h_y * j < h)
                T[i * M + j] = T_1;
            else
                T[i * M + j] = T_2;
        }
}
int main() {
    int a = 1, i, j, time = 5000, t;
    float ai, bi, ci, fi;
    float h_x, h_y, tau, n;
    double **alpha = new double * [2];
    // for mpi implementation:
    int size, rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *T = new double[M * M];
    double *T_proc = new double[M * M/size];

    ofstream on("file_mpi.txt");
    /*on << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"y\", \"x\", \"T\"" << endl <<
    "ZONE T = \"Numerical\", I = " << N_x << ", J = " << N_y << ", F = Point";*/

    if(!on){
        cout << "Error openning input file. \n";
        return -1;
    }

    for (i = 0; i < 2; i++)
        alpha[i] = new double [2];
    /*for (i = 0; i < N_x; i++){
        T[i] = new double *[N_y/size];*/
       /* for (j = 0; j < N_y; j++){
            T[i][j] = new double [time];
        }
    }*/

    h_x = L_x/(N_x - 1);
    h_y = L_y/(N_y - 1);
    tau = 1/(pow(h_x,-2)+pow(h_y,-2))*0.5; // ?????????????????????? ??????????
    // ???????????????? ????????????????????????: tau <= h^2/(2*p), p - ?????????? ??????
    // ????. ??????????????????, ?????????? "???????????????????????? ???????????????????? ????????" ??????. 314
    initMat(T, h_x, h_y);

    for(i = 0; i < M; i++) {
        for (j = 0; j < M; j++){
            on << T[i * M + j] << ' ';
        }
        on << endl;
    }

    /*for (t = 0; t < time - 1; t++){
        for (i = 0; i < N_x; i++){
            for (j = 0; j < N_y; j++) {
                if (j > 0 && j < N_y-1 && i == 0) {
                    alpha[0][1] = 0;
                    alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                    alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                    alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                }
                else if (j == N_y-1 && i == 0) {
                    alpha[0][1] = 0;
                    alpha[1][0] = 0;
                    alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                    alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                }
                else if (j > 0 && j < N_y-1 && i == N_x-1){
                    alpha[0][0] = 0;
                    alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                    alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                    alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                }
                else if (i > 0 && i < N_x-1 && j == 0) {
                    alpha[1][1] = 0;
                    alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                    alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                    alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                }
                else if (i == N_x-1 && j == 0) {
                    alpha[0][0] = 0;
                    alpha[1][1] = 0;
                    alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                    alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                }
                else if (i > 0 && i < N_x-1 && j == N_y-1) {
                    alpha[1][0] = 0;
                    alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                    alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                    alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                }
                else if (i == 0 && j == 0) {
                    alpha[0][1] = 0;
                    alpha[1][1] = 0;
                    alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                    alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                }
                else if (i == N_x-1 && j == N_y-1) {
                    alpha[0][0] = 0;
                    alpha[1][0] = 0;
                    alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                    alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                }
                else {
                    alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                    alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                    alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                    alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                }
                T[i][j][t+1] = tau * ((alpha[0][0] + alpha[0][1])/pow(h_x,2) + (alpha[1][0] + alpha[1][1])/pow(h_y,2)) + T[i][j][t];
            }
        }
    }*/
    cout << "ok" << endl;

    /*for (j = 0; j < N_y; j++){
        for(i = 0; i < N_x; i++){
            on << endl;
            on << i << " " << j << " " << T[i][j][4500];
        }*/
        //on << endl;
    //}

    delete[] T;
    delete[] alpha;

    MPI_Finalize();
    return 0;
}