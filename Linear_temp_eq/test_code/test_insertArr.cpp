#include <iostream>
#include <math.h>
#define PI 3.14

using namespace std;

void InsertArr(double x0[2][6], int ind, int N_total, int N_shift, char x){
    int i, j;
    for (i = N_shift; i >= ind; i--)
        for (j = 0; j < N_total; j++){
            if (x == 'x')
                x0[i][j] = x0[i-1][j];
            else if (x == 'y')
                x0[j][i] = x0[j][i-1];
        }
}


int main(){
    double mass[2][6] = {{1.1, 4.3, 6.7, 7.0, 0.0, 0.0},{3.7, 8.9, 3.4, 9.0, 0.0, 0.0}};
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 5; j++)
            cout << mass[i][j] << " ";
        cout << endl;
    }

    InsertArr(mass, 2, 2, 4, 'y');

    InsertArr(mass, 3, 2, 5, 'y');

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 5; j++)
            cout << mass[i][j] << " ";
        cout << endl;
    }
    /*for (i = N_x/x_domains; i < N_x; i+=N_x/x_domains){
        InsertArr(x0, i+cnt, N_y_total, cnt + N_x_global, 'x');
        cnt += 1;
        InsertArr(x0, i + 2 + cnt, N_y_total, cnt + N_x_global, 'x');
        cnt += 1;
    }
    cnt = 0;
    for (j = N_y/y_domains; j < N_y; j+=N_y/y_domains){
        InsertArr(x0, j + cnt, N_x_total, cnt + N_y_global, 'y');
        cnt += 1;
        InsertArr(x0, j + 2 + cnt, N_x_total, cnt + N_y_global, 'y');
        cnt += 1;
    }
    InsertArr(x0, 1, N_y_total, N_x_total - 2, 'x');
    InsertArr(x0, N_x_total - 2, N_y_total, N_x_total - 1, 'x');

    InsertArr(x0, 1, N_x_total, N_y_total - 2, 'y');
    InsertArr(x0, N_y_total - 2, N_x_total, N_y_total - 1, 'y');*/
    return 0;
}