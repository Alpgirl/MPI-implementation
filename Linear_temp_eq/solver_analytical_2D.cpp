#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#define PI 3.14
using namespace std;

int main() {
    int i, j, N_x = 4, N_y = 4, time = 1, t, t_c = 0, z;
    cout << setprecision(15);

    double L_x = 10.0, L_y = 20.0, T_1 = 1000.0, T_2 = 300.0, ai, bi, ci, fi;
    double h_x, h_y, tau, n, h = 7.0, coef_time, a_2 = 0.5;
    double test, test1, test2;
    double ***T = new double ** [(N_x)];
    double **alpha = new double * [2];
    double **T_exact = new double * [N_x];

    ofstream on_exact_1("file_analytical_compare_mpi_exact.dat");
    on_exact_1 << "TITLE = \"Analytic sol\"" << endl << "VARIABLES = \"x\", \"T\"" << endl <<
    "ZONE T = \"Numerical\", I = " << N_x << ", F = Point" << endl;

    ofstream on_exact_2("file_analytical_compare_mpi_calc.dat");
    /*on_exact_2 << "TITLE = \"Analytic sol\"" << endl << "VARIABLES = \"x\", \"T\"" << endl <<
    "ZONE T = \"Numerical\", I = " << N_x << ", F = Point" << endl;*/


    if(!on_exact_1 || !on_exact_2){
        cout << "Error openning input file. \n";
        return -1;
    }

    for (i = 0; i < 2; i++)
        alpha[i] = new double [2];
    for (i = 0; i < N_x; i++){
        T[i] = new double *[N_y];
        T_exact[i] = new double[N_y];
        for (j = 0; j < N_y; j++){
            T[i][j] = new double [time];
        }
    }

    h_x = L_x/(N_x-1);
    h_y = L_y/(N_y-1);
    tau = 1.0/4*pow(min(h_x,h_y),2)/a_2; // оптимальное время
    // критерий устойчивости: tau <= h^2/(2*p), p - число мер
    // см. Самарский, Гулин "Устойчивость разностных схем" стр. 314

    // Задаем начальное температурное поле
    // Распределение по прямоугольникам

    /*for(i = 0; i < N_x; i++)
        for (j = 0; j < N_y; j++){
            if (i < h/h_x - 1 and j < h/h_y - 1)
                T[i][j][0] = T_1;
            else
                T[i][j][0] = T_2;
        }*/
    

    // Начальные условия для проверки задачи
    for (i = 0; i < N_x; i++)
        for(j = 0; j < N_y; j++){
            if(i == 0 or i == N_x - 1 or j == 0 or j == N_y - 1) T[i][j][0] = 0;
            else {T[i][j][0] = sin(PI*(h_x*i))*sin(PI*(h_y*j)); /*cout << h_x << " " << i << " " << j << endl;*/}
        }
    // Аналитическое решение
    coef_time = exp(-pow(PI,2)*t_c*tau);
    for (i = 0; i < N_x; i++){
        for (j = 0; j < N_y; j++){
            if(i == 0 or i == N_x - 1 or j == 0 or j == N_y - 1) T_exact[i][j] = 0;
            else T_exact[i][j] = coef_time * sin(PI*(h_x*i))*sin(PI*(h_y*j));
        }
    }

    for (t = 0; t < time - 1; t++){
        for (i = 0; i < N_x; i++){
            for (j = 0; j < N_y; j++) {
                if(i == 0 or i == N_x-1 or j == 0 or j == N_y - 1) T[i][j][t] = 0;
                else {
                    if (j > 0 && j < N_y-1 && i == 0) {
                        alpha[0][1] = 0.0;
                        alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                        alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                        alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                    }
                    else if (j == N_y-1 && i == 0) {
                        alpha[0][1] = 0.0;
                        alpha[1][0] = 0.0;
                        alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                        alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                    }
                    else if (j > 0 && j < N_y-1 && i == N_x-1){
                        alpha[0][0] = 0.0;
                        alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                        alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                        alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                    }
                    else if (i > 0 && i < N_x-1 && j == 0) {
                        alpha[1][1] = 0.0;
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
                        alpha[1][0] = 0.0;
                        alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                        alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                        alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                    }
                    else if (i == 0 && j == 0) {
                        alpha[0][1] = 0.0;
                        alpha[1][1] = 0.0;
                        alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                        alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                    }
                    else if (i == N_x-1 && j == N_y-1) {
                        alpha[0][0] = 0.0;
                        alpha[1][0] = 0.0;
                        alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                        alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                    }
                    else {
                        alpha[0][0] = T[i+1][j][t] - T[i][j][t];
                        alpha[0][1] = - T[i][j][t] + T[i-1][j][t];
                        alpha[1][0] = T[i][j+1][t] - T[i][j][t];
                        alpha[1][1] = - T[i][j][t] + T[i][j-1][t];
                    }
                    T[i][j][t+1] = tau * a_2 * ((alpha[0][0] + alpha[0][1])/pow(h_x,2) + (alpha[1][0] + alpha[1][1])/pow(h_y,2)) + T[i][j][t];
                }
            }
        }
    }
    cout << "ok" << endl;
    cout << h_x << " " << h_y << " " << tau << endl;

    for (j = 0; j < N_y; j++){
        for (i = 0; i < N_x; i++){
            //on_exact_1 << endl;
            //on_exact_1 << i << " " << T_exact[i][10] << endl;// << " ";  /*T[i][j][t_c] << " ";*/
            on_exact_2 <</* i << " " <<*/ T[i][j][t_c] << " ";
        }
        //on_exact_1 << endl;
        on_exact_2 << endl;
    }
    delete[] T;
    delete[] alpha;
    return 0;

}
