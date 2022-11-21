#include <iostream>
#include <fstream>
#include <cmath>
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

int main() {
    int a = 1, i, j, N_x = 2048, N_y = 2048, time = 51, t;
    float L_x = 10.0, L_y = 20.0, T_1 = 1000.0, T_2 = 300.0, ai, bi, ci, fi;
    float h_x, h_y, tau, n, h = 5.0;
    double ***T = new double ** [(N_x)];
    double **alpha = new double * [2];

    ofstream on("file_compare_with_mpi.dat");
    /*on << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"y\", \"x\", \"T\"" << endl <<
    "ZONE T = \"Numerical\", I = " << N_x << ", J = " << N_y << ", F = Point";*/

    if(!on){
        cout << "Error openning input file. \n";
        return -1;
    }

    for (i = 0; i < 2; i++)
        alpha[i] = new double [2];
    for (i = 0; i < N_x; i++){
        T[i] = new double *[N_y];
        for (j = 0; j < N_y; j++){
            T[i][j] = new double [time];
        }
    }

    h_x = L_x/(N_x+2);
    h_y = L_y/(N_y+2);
    tau = 1./(4*pow(a,2))*pow(min(h_x,h_y),2); // оптимальное время
    // критерий устойчивости: tau <= h^2/(2*p), p - число мер
    // см. Самарский, Гулин "Устойчивость разностных схем" стр. 314

    for(i = 0; i < N_x; i++)
        for (j = 0; j < N_y; j++){
            if (i < h/h_x - 1 and j < h/h_y - 1)
                T[i][j][0] = T_1;
            else
                T[i][j][0] = T_2;
        }

    /*for (j = 0; j < N_y - 1; j++){
        for(i = 0; i < N_x - 1; i++){
            on << T[i][j] << ' ';
        }
        on << endl;
    }*/

    for (t = 0; t < time - 1; t++){
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
                T[i][j][t+1] = pow(a,2) * tau * ((alpha[0][0] + alpha[0][1])/pow(h_x,2) + (alpha[1][0] + alpha[1][1])/pow(h_y,2)) + T[i][j][t];
            }
        }
    }
    cout << "ok" << endl;
    //on << endl;
    for (j = 0; j < N_y; j++){
        for(i = 0; i < N_x; i++){
            //on << endl;
            on << T[i][j][50] << " ";
        }
        on << endl;
    }

    delete[] T;
    delete[] alpha;
    return 0;

}