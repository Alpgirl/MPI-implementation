#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <cstdlib>
#include "updateBound.h"
#define h  5.0
#define L_x  10.0
#define L_y  20.0
using namespace std;

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
void initValues(double **x0, int N_x_total, int N_y_total, double T_1, double T_2, 
                double h_x, double h_y, double x_domains, double y_domains, int N_x_global, int N_y_global, int N_x, int N_y){
    int i, j, cnt;
    double T_temp;
    double *x_index = new double[N_x_total];
    double *y_index = new double[N_y_total];

    for (i = 0; i < N_x_global; i++)
        for (j = 0; j < N_y_global; j++){
            if (i < h/h_x and j < h/h_y) x0[i][j] = T_1;
            else x0[i][j] = T_2;
        }
    cnt = 0;
    for (i = N_x/x_domains; i < N_x; i+=N_x/x_domains){
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
    InsertArr(x0, N_y_total - 2, N_x_total, N_y_total - 1, 'y');
}
void ProcessToMap(int *xs, int *ys, int *xe, int *ye, int xcell, int ycell, int x_domains, int y_domains){
    int i, j;
    // computation of starting ys, ye on (Ox) standard axis for the first column of global domain
    // convention (i,j) = (row,column)
    for (i = 0; i < x_domains; i++){
        ys[i] = 2;
        ye[i] = ys[i] + ycell - 1;
    }
    // computation of ys, ye on (Ox) standard axis for other columns of global domain
    for (i = 1; i < y_domains; i++){
        for (j = 0; j < x_domains; j++){
            ys[i * x_domains + j] = ys[(i - 1) * x_domains + j] + ycell + 2;
            ye[i * x_domains + j] = ys[i * x_domains + j] + ycell - 1;
        }
    }
    // computation of starting xs, xe on (Oy) standard axis for the first column of global domain
    for (i = 0; i < y_domains; i++){
        xs[i * x_domains] = 2;
        xe[i * x_domains] = xs[i * x_domains] + xcell -1;
    }
    // computation of xs, xe on (Oy) standard axis for other columns of global domain
    for (i = 1; i <= y_domains; i++){
        for (j = 1; j < x_domains; j++){
            xs[(i - 1) * x_domains + j] = xs[(i - 1)* x_domains + (j - 1)] + xcell + 2;
            xe[(i - 1) * x_domains + j] = xs[(i - 1) * x_domains + j] + xcell - 1;
        }
    }
}
int main() {
    int i, j, k, l, n, t, zx = 0, zy = 0, time = 400;

    // diffusivity coefficient
    double a = 1.0;

    // physical parameters
    double T_1 = 1000.0, T_2 = 300.0;

    // time variables
    double time_init, time_final, elapsed_time;

    // Various parameters for dimensions
    int N_x = 102, N_y = 102, x_domains = 2, y_domains = 2;
    int N_x_global, N_y_global;
    int N_x_total, N_y_total;

    // space and time steps
    double tau, h_x, h_y;

    // for mpi implementation:
    int size, rank, ndims;
    MPI_Comm comm, comm2d;
    int dims[2];
    int periods[2];
    int reorganisation = 0;
    int S = 0, E = 1, N = 2, W = 3;
    int neighBor[4];
    int xcell, ycell;

    // MPI initialisation
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ((size == 0) or size != (x_domains * y_domains)) {
        cout << "Number of processes not equal to Number of subdomains" << endl;
    }

    N_x_global = N_x + 2;
    N_y_global = N_y + 2;

    h_x = L_x/N_x_global;
    h_y = L_y/N_y_global;
    tau = 1./(4.0*pow(a,2))*pow(min(h_x,h_y),2); // оптимальное время
    // критерий устойчивости: tau <= h^2/(2*p), p - число мер
    // см. Самарский, Гулин "Устойчивость разностных схем" стр. 314
    // tau = 1/4a^2 * min(h_x,h_y)^2

    N_x_total = N_x + 2 * x_domains + 2;
    N_y_total = N_y + 2 * y_domains + 2;

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

    // create 2D cartesion grid
    periods[0] = 0;
    periods[1] = 0;
    
    // number of dimensions
    ndims = 2;

    // invert (Ox, Oy) classic convention
    dims[0] = y_domains;
    dims[1] = x_domains;

    // makes a new communicator to which topology information has benn attached
    MPI_Cart_create(comm, ndims, dims, periods, reorganisation, &comm2d);

    // indentify neighBors; specify a "dummy" source or destination for communication
    neighBor[0] = MPI_PROC_NULL;
    neighBor[1] = MPI_PROC_NULL;
    neighBor[2] = MPI_PROC_NULL;
    neighBor[3] = MPI_PROC_NULL;

    // "фиктивный" отправитель/получатель. Необходимо для работы с границами.

    // left/west and right/east neighBors; returns the shifted source and destination ranks, given a shft direction and amount
    MPI_Cart_shift(comm2d, 0, 1, &neighBor[W], &neighBor[E]);

    // bottom/sourth and upper/north neighBors; -||-
    MPI_Cart_shift(comm2d, 1, 1, &neighBor[S], &neighBor[N]);

    // size of each cell
    xcell = N_x/x_domains;
    ycell = N_y/y_domains;

    double *xtemp = new double [xcell * ycell];

    // compute xs, ys, xe, ye for each cell on the grid
    ProcessToMap(xs, ys, xe, ye, xcell, ycell, x_domains, y_domains);

    // initialize values
    initValues(x0, N_x_total, N_y_total, T_1, T_2, h_x, h_y, x_domains, y_domains, N_x_global, N_y_global, N_x, N_y);

    updateBound(x0, neighBor, comm2d, rank, xs, ys, xe, ye, ycell);

    /* Starting time */
    time_init = MPI_Wtime();
    
    for (t = 0; t < time; t++){
        // solver
        computeNext(x0, x, tau, h_x, h_y, rank, xs, ys, xe, ye, a);

        // communication between processes
        updateBound(x0, neighBor, comm2d, rank, xs, ys, xe, ye, ycell);

        // update external boundaries of subdomains
        for (n = 0; n < x_domains; n++){
            for (i = xs[n]-1; i <= xe[n]+1; i++)
                x0[i][ys[n]-1] = x0[i][ys[n]];
            for (i = xs[size - n - 1]-1; i <= xe[size - n - 1]+1; i++)
                x0[i][ye[size - n - 1]+1] = x0[i][ye[size - n - 1]];
        }
        for (n = 1; n <= y_domains; n++){
            for (j = ys[x_domains * n - 1]-1; j <= ye[x_domains * n - 1]+1; j++)
                x0[xe[x_domains * n - 1]+1][j] = x0[xe[x_domains * n - 1]][j];
            for (j = ys[x_domains * (n - 1)]-1; j <= ye[x_domains * (n - 1)]+1; j++)
                x0[xs[x_domains * (n - 1)]-1][j] = x0[xs[x_domains * (n - 1)]][j];
        }
        // printing for testing working process of each process
        if (rank == 0 and t == time - 1) {
            ofstream test("just_test_mpi_0.dat");
            for (i = 0; i < N_x_total; i++){
                for (j = 0; j < N_y_total; j++){ 
                    test << x0[i][j] << " ";
                }
                test << endl;
            }
            test.close();
        }
        if (rank == 1 and t == time - 1) {
            ofstream test("just_test_mpi.dat");
            for (i = 0; i < N_x_total; i++){
                for (j = 0; j < N_y_total; j++){ 
                    test << x0[i][j] << " ";
                }
                test << endl;
            }
            test.close();
        }
    }

    // Ending time 
    time_final = MPI_Wtime();

    // Elapsed time
    elapsed_time = time_final - time_init;

    j=1;
    for (i=xs[rank];i<=xe[rank];i++) {
        for (k=0;k<ycell;k++){
            xtemp[(j-1)*ycell+k] = x0[i][ys[rank]+k];
        }
        j=j+1;
    }

    MPI_Gather(xtemp, xcell*ycell , MPI_DOUBLE , xfinal, xcell*ycell, MPI_DOUBLE, 0 , comm);

    //cout << ys[rank] << " " << ye[rank] << endl;
    //cout << xs[rank] << " " << xe[rank] << endl;
    //cout << N_x << " " << N_x_total << " " << N_x_global << endl;

    if (rank == 0){
        ofstream on("file_mpi_vizual.dat");
        on << "TITLE = \"Bivariate normal distribution density\"" << endl << "VARIABLES = \"y\", \"x\", \"T\"" << endl <<
        "ZONE T = \"Numerical\", I = " << N_x << ", J = " << N_y << ", F = Point" << endl;
        ofstream on2("file_mpi.dat");

        if(!on or !on2){
            cout << "Error openning input file. \n";
            return -1;
        }

        for (i=1;i<=y_domains;i++){
            for (j=0;j<ycell;j++) {
                for (k=1;k<=x_domains;k++) {
                    for (l=0;l<xcell;l++){
                        on << zx << " " << zy << " " << xfinal[(i-1)*x_domains*xcell*ycell+(k-1)*xcell*ycell+l*ycell+j] << endl;
                        on2 << xfinal[(i-1)*x_domains*xcell*ycell+(k-1)*xcell*ycell+l*ycell+j] << " ";
                        zy++;
                    }
                }
                zy = 0;
                on2 << endl;
                zx++;
            }
        }
        on.close();
        on2.close();
        cout << endl;
        cout << tau << endl;
        cout << elapsed_time << endl;
    }

    delete[] x0[0];
    delete[] x0;
    delete[] x;
    delete[] x[0];
    delete[] xtemp;
    delete[] xfinal;
    delete[] xs;
    delete[] ys;
    delete[] xe;
    delete[] ye;
    MPI_Finalize();
    return 0;
}