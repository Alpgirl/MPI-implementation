#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include <mpi.h>
#include "updateBound_stick.h"

#define h  7.0
#define L_x  0.5//10.0
#define L_y  0.5//20.0
#define ro   8960. // density
#define lbd  401. // coefficient of thermal conductivity
#define c    380. // heat capacity
#define PI 3.14

using namespace std;

void fill_null_1d(int *x0, int N){
    for (int i = 0; i < N; i++)
        x0[i] = 0.0;
}



void fill_null_2d(double **x0, int N1, int N2){
    for (int i = 0; i < N1; i++)
        for(int j = 0; j < N2; j++)
            x0[i][j] = 0.0;
}



void InsertArr(double **x0, int ind, int N_total, int N_shift, char x){
    int i, j;
    for (i = N_shift; i >= ind; i--)
        for (j = 0; j < N_total; j++){
            if (x == 'x')
                x0[i][j] = x0[i-1][j];
            else if (x == 'y')
                x0[j][i] = x0[j][i-1];
        }
}


void PrepareParalField (double **x0,  int N_x_total, int N_y_total,
                    int N_x_global, int N_y_global, int N_x, int N_y, int size, int* dims, int rank){
    int i, j = 0, cnt;
    cnt = 0; 
    for (i = 1; i < size; i++){
        j += dims[i];
        InsertArr(x0, j+cnt+1, N_x_total, cnt + N_y_global, 'y');
        cnt += 1;
        InsertArr(x0, j + 2 + cnt, N_x_total, cnt + N_y_global, 'y');
        cnt += 1;
    }
    InsertArr(x0, 1, N_y_total, N_x_total - 2, 'x');
    InsertArr(x0, N_x_total - 1, N_y_total, N_x_total - 1, 'x');

    InsertArr(x0, 1, N_x_total, N_y_total - 2, 'y');
    InsertArr(x0, N_y_total - 1, N_x_total, N_y_total - 1, 'y');
}



void initValueStick(double **x0, int N_x_total, int N_y_total,
                    int N_x_global, int N_y_global, int N_x, int N_y, int size, int* dims, int rank,
                    int T_l = 80, int T_r = 30, int T_u = 100, int T_d = 0, int T_i = 5){
    int i, j;

    for (i = 0; i < N_x_global; i++)
        for (j = 0; j < N_y_global; j++){
            if(i == 0) x0[i][j] = T_l;
            else if(i == N_x_global-1) x0[i][j] = T_r;
            else if(j == 0) x0[i][j] = T_u;
            else if(j == N_y_global-1) x0[i][j] = T_d;
            else x0[i][j] = T_i;
        }
    PrepareParalField(x0, N_x_total, N_y_total, N_x_global, N_y_global, N_x, N_y, size, dims, rank);
}


void ProcessToMap(int *xs, int *ys, int *xe, int *ye, int *dims, int size, int N_x){
    int i, j, proc;
    // computation of starting ys, ye on (Ox) standard axis for the first column of global domain
    // convention (i,j) = (row,column)
    for (proc = 0; proc < size; proc++){
        xs[proc] = 2;
        xe[proc] = xs[proc] + N_x - 1;
    }
    ys[0] = 2;
    ye[0] = ys[0] + dims[0] - 1;
    for (proc = 1; proc < size; proc++){
        ys[proc] = ye[proc - 1] + 3;
        ye[proc] = ys[proc] + dims[proc] - 1;
    }
}



int main() {
    int i, j, k, l, n, t, cnt = 0, time_step = 10000, time_to = 60;
    double zx = 0, zy = 0, T;
    cout << setprecision(15);

    // diffusivity coefficient
    double a_2 = lbd/(ro * c) ;//0.5;

    // time variables
    double time_init, time_final, elapsed_time;

    // Various parameters for dimensions
    int N_x = 150, N_y = 150, Ncell, remainCell;
    int N_x_global, N_y_global;
    int N_x_total, N_y_total;

    // space and time steps
    double tau, h_x, h_y;

    // for mpi implementation:
    int size, rank, ndims;
    MPI_Comm comm, comm2d;

    // MPI initialisation
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    
    int* dims = new int[size];
    int* dims_shift = new int[size];
    int* dims_recv = new int[size];

    if (size == 0) {
        cout << "Number of processes not equal to Number of subdomains" << endl;
    }

    // Calculate number of cells for each process and save information in dims
    remainCell = N_y - N_y/size * size;

    for (i = 0; i < size; i++){
        dims[i] = N_y/size;
    }

    for (i = 0; i < remainCell; i++){
        if (cnt == size - 1) cnt = 0;
        dims[cnt] += 1;
        cnt++;
    }
    dims_shift[0] = 0;
    dims_recv[0] = N_x * dims[0];
    for (i = 1; i < size; i++){
        dims_shift[i] = dims[i-1] * N_x + dims_shift[i-1];
        dims_recv[i] = N_x * dims[i];
    }

    N_x_global = N_x + 2;
    N_y_global = N_y + 2;

    h_x = L_x/(N_x_global-3);
    h_y = L_y/(N_y_global-3);
    tau = 1./4.0*pow(min(h_x,h_y),2)/a_2; // оптимальное время
    // критерий устойчивости: tau <= h^2/(2*p), p - число мер
    // см. Самарский, Гулин "Устойчивость разностных схем" стр. 314
    // tau = 1/4a^2 * min(h_x,h_y)^2

    N_x_total = N_x_global + 2;
    N_y_total = N_y_global + 2 * size;

    // Allocate 2D contiguous arrays x and x0 */
    double **x0 = new double * [N_x_total];
    double **x = new double * [N_x_total];
    double *xfinal = new double [N_x * N_y];

    for (i = 0; i < N_x_total; i++){
        x0[i] = new double[N_y_total];
        x[i] = new double[N_y_total];
    }

    // allocate coordinates of processes
    int *xs = new int[size];
    int *xe = new int[size];
    int *ys = new int[size];
    int *ye = new int[size];

    fill_null_1d(xs, size);
    fill_null_1d(xe, size);
    fill_null_1d(ys, size);
    fill_null_1d(ye, size);
    fill_null_2d(x0, N_x_total, N_y_total);
    fill_null_2d(x, N_x_total, N_y_total);

    double *xtemp = new double [N_x * dims[rank]];

    // compute xs, ys, xe, ye for each cell on the grid
    ProcessToMap(xs, ys, xe, ye, dims, size, N_x);

    // initialize values

    initValueStick(x0, N_x_total, N_y_total, N_x_global, N_y_global, N_x, N_y, size, dims, rank);


    /* Starting time */
    time_init = MPI_Wtime();
    
    for (t = 0; t < time_step; t++){
        if(tau * t > time_to) break;
        // solver
        computeNext(x0, x, tau, h_x, h_y, rank, xs, ys, xe, ye, a_2);

        // communication between processes
        updateBound2(x0, N_x_global, comm, rank, size, xs, ys, xe, ye);
    }
    
    // Ending time 
    time_final = MPI_Wtime();

    // Elapsed time
    elapsed_time = time_final - time_init;

    j=0;
    for (k=ys[rank];k<=ye[rank];k++){
        for (i=xs[rank];i<=xe[rank];i++) {
            xtemp[j] = x0[i][k];
            j++;
        }
    }


    MPI_Gatherv(&xtemp[0], N_x*dims[rank] , MPI_DOUBLE, &xfinal[0], dims_recv , dims_shift, MPI_DOUBLE , 0, comm);

    if (rank == 0){
        ofstream on("file_mpi_vizual.dat");
        on << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"x\", \"y\", \"T\"" << endl <<
        "ZONE T = \"Numerical\", I = " << N_x << ", J = " << N_y << ", F = Point" << endl;
        ofstream on2("file_mpi_analyt.dat");
        ofstream on1d_016x("TECPLOT for report/check 1D solution/file_num_y = 0.16_t=60.dat");
        on1d_016x << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"x\", \"T\"" << endl <<
        "ZONE T = \"Numerical\", J = " << N_x + 2<< ", F = Point" << endl;

        ofstream on1d_016y("TECPLOT for report/check 1D solution/file_num_x = 0.16_t=60.dat");
        on1d_016y << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"y\", \"T\"" << endl <<
        "ZONE T = \"Numerical\", J = " << N_y + 2<< ", F = Point" << endl;

        ofstream on1d_032x("TECPLOT for report/check 1D solution/file_num_y = 0.32_t=60.dat");
        on1d_032x << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"x\", \"T\"" << endl <<
        "ZONE T = \"Numerical\", J = " << N_x + 2<< ", F = Point" << endl;

        ofstream on1d_032y("TECPLOT for report/check 1D solution/file_num_x = 0.32_t=60.dat");
        on1d_032y << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"y\", \"T\"" << endl <<
        "ZONE T = \"Numerical\", J = " << N_y + 2<< ", F = Point" << endl;

        if(!on or !on2){
            cout << "Error openning input file. \n";
            return -1;
        }

        on1d_016x << 0 << " " << 80.0 << endl;
        on1d_032x << 0 << " " << 80.0 << endl;
        on1d_016y << 0 << " " << 100.0 << endl;
        on1d_032y << 0 << " " << 100.0 << endl;
        

        for (k = 0; k < N_y; k++){
            for (j = 0; j < N_x; j++){
                zx = h_x * (j+1);
                zy = h_y * (k+1);
                if (j == round(0.16/h_x)) {
                    on1d_016y << zy << " " << xfinal[j + k*N_x] << endl;
                }
                if (k == round(0.16/h_y)) {
                    on1d_016x << zx << " " << xfinal[j + k*N_x]<< endl;
                }
                if (j == round(0.32/h_x)) {
                    on1d_032y << zy << " " << xfinal[j + k*N_x] << endl;
                }
                if (k == round(0.32/h_y)) {
                    on1d_032x << zx<< " " << xfinal[j + k*N_x]  << endl;
                }
                on2 << xfinal[j + k*N_x] << " ";
                on << h_x * j << " " << h_y * (N_y - k - 1) << " " << xfinal[j + k*N_x] << endl;
            }
            on2 << endl;
        }

        on1d_016x << (N_x+1) * h_x << " " << 30.0 << endl;
        on1d_032x << (N_x+1) * h_x << " " << 30.0 << endl;
        on1d_016y << (N_y+1 ) * h_y << " " << 0.0 << endl;
        on1d_032y << (N_y+1) * h_y << " " << 0.0 << endl;
        

        on.close();
        on2.close();
        cout << endl;
        cout << "Time steps: " << t << endl;
        cout << "Elapsed time: " << elapsed_time << endl;
        cout << "x step: " << h_x << ", y step: " << h_y << ", time step: " << tau << endl;
    }
    for (i = 0; i < N_x_total; i++){
        delete x0[i];
        delete x[i];
    }
    delete[] x0;
    delete[] x;
    delete[] xtemp;
    delete[] xfinal;
    delete[] xs;
    delete[] ys;
    delete[] xe;
    delete[] ye;
    MPI_Finalize();
    return 0;
}