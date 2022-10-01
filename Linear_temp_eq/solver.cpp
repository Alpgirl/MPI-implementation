#include "solver.h"
#include <iostream>
using namespace std;

int main() {
    int a = 1, i, j;
    float L_x = 10.0, L_y = 20.0, T_1 = 300.0, T_2 = 270.0, N_x = 100, N_y = 100, time = 100, ai, bi, ci, fi, alpha, betta;
    float h_x, h_y, tau, n, h = 5.0;
    double **T = new double* [N_x + 2]; 

    for (i = 0; i < N_x + 2; i++){
        T[i] = new double [N_y + 2];
    }

    h_x = L_x/(N_x - 1);
    h_y = L_y/(N_y - 1);
    tau = time/100.0;

    for(i = 1, i < N_x, i++)
        for (j = 1, j < N_y, j++){
            if (h_x * i <= h && h_y * j <= h)
                T[i][j] = T1;
            else
                T[i][j] = T2;
        }

    delete T;

}