// USEFUL REFERENCES///
// https://www.boost.org/doc/libs/master/doc/html/mpi/tutorial.html //

#include <iostream>
#include <vector>
#include <mpi.h>

int main(int argc, char* argv[]) {
	int rank, size, tag = 1;
	double c = 1.;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	// Send-Recv with 1d data type
	/*if (rank == 0) {
		std::vector <double> x = { 0.1, 3.5, -7., 8.9 };
		int n = x.size();
		if (rank == 0) std::cout << "Proc with rank " << rank << " says: " << x.size();

		MPI_Send(&n, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
		MPI_Send(x.data(), n, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
	}

	else {
		int n;
		MPI_Recv(&n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		std::vector<double> y(n);
		MPI_Recv(y.data(), n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

		std::cout << "Processor " << rank << " received" << n << " pieces of data: ";
		for (double f : y) std::cout << f << " ";
		std::cout << '\n';
	}*/

	// Send-Recv with 2d data type
	if (rank == 0) {
		std::vector < std::vector <double> > x;//(3, std::vector <double>(5));
		for (int i = 0; i < 3; i++) {
			std::vector <double> vv;
			for (int j = 0; j < 5; j++) {
				vv.push_back(j * c);
				c += 0.1;
			}
			x.push_back(vv);
		} 
		for (int i = 0; i < x.size(); i++) {
			for (int j = 0; j < x[i].size(); j++)
				std::cout << x[i][j] << " ";
			std::cout << '\n';
		}

		int n1 = x[1].size();  std::cout << "Proc " << rank << " says col: " << n1 << std::endl;
		int n2 = x.size();	 std::cout << "Proc " << rank << " says row: " << n2 << std::endl;

		MPI_Send(&n1, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
		MPI_Send(&n2, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
		for (int i = 0; i < 3; i++)
			MPI_Send(x[i].data(), n1, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
	}

	else {
		int n1, n2;
		MPI_Recv(&n1, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&n2, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		std::vector < std::vector <double> > y(n2, std::vector <double> (n1));
		for (int i = 0; i < n2; i++)
			MPI_Recv(y[i].data(), n1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

		std::cout << "Processor " << rank << " received " << n1*n2 << " pieces of data: \n";
		for (int i = 0; i < y.size(); i++) {
			for (int j = 0; j < y[i].size(); j++)
				std::cout << y[i][j] << " ";
			std::cout << '\n';
		}
	}

	MPI_Finalize();
	return 0;
}