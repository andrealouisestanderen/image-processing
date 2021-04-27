#include <assert.h>
#include <cmath>
#include "../lab3/png++-0.2.7/png.hpp"
#include "stdio.h"
#include "string.h"
#include <string>
#include <sstream>
#include <chrono>
#include <omp.h>
#include "mpi.h"

using namespace std;

typedef vector<double> Array;
typedef vector<Array> Matrix;
typedef vector<Matrix> Image;

int main(int argc, char *argv[])
{

    auto t1 = std::chrono::high_resolution_clock::now();

    auto t1_1 = std::chrono::high_resolution_clock::now();

    int rank, size, tag, rc, i, j, row_id;
    double *send, *recv, *send_back, *recv_back, data;
    MPI_Status status;
    Array row, newRow;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Matrix testMatrix = {{1, 1, 1},
                         {2, 2, 2},
                         {3, 3, 3}};

    int row_size = testMatrix[0].size();
    int col_size = testMatrix.size();
    int batch_size = row_size / (size - 1);

    send = (double *)malloc(sizeof(double) * row_size);
    recv = (double *)malloc(sizeof(double) * row_size);
    send_back = (double *)malloc(sizeof(double) * row_size);
    recv_back = (double *)malloc(sizeof(double) * row_size);

    if (rank == 0)
    {
        int row_start = 0;
        for (i = 1; i < size; i++)
        {
            for (int row_id = row_start; row_id < i * batch_size; row_id++)
            {
                for (int col_id = 0; col_id < col_size; col_id++)
                {
                    double data = testMatrix[row_id][col_id];
                    send[col_id] = data;
                    std::cout << "\nAllocating data into mem:\n"
                              << send[col_id];
                }
                rc = MPI_Send(send, row_size, MPI_DOUBLE, i, row_id, MPI_COMM_WORLD);
            }
            row_start = i * batch_size;
        }
    }
    else
    {

        for (int row_id = (rank-1) * batch_size; row_id < (rank-1) * batch_size + batch_size; row_id++)
        {
            rc = MPI_Recv(recv, row_size, MPI_DOUBLE, 0, row_id, MPI_COMM_WORLD, &status);
            for (int col_id = 0; col_id < col_size; col_id++)
            {
                double data = recv[col_id] * 2;
                std::cout << "\nAltering received data: \n"
                          << data;
                send_back[col_id] = data;
            }
            rc = MPI_Send(send_back, row_size, MPI_DOUBLE, 0, row_id, MPI_COMM_WORLD);
        }
    }

    if (rank == 0)
    {
        for (i = 1; i < size; i++)
        {
            for (int row_id = i * batch_size - 1; row_id < i * batch_size + batch_size - 1; row_id++)
            {
                rc = MPI_Recv(recv_back, row_size, MPI_DOUBLE, i, row_id, MPI_COMM_WORLD, &status);
                for (int col_id = 0; col_id < col_size; col_id++)
                {
                    double data = recv_back[col_id];
                    std::cout << "\nDATA RECEIVED BACK AT MASTER: \n"
                              << data;
                }
            }
        }
    }

    auto t2_1 = std::chrono::high_resolution_clock::now();
    auto duration_1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2_1 - t1_1).count();
    std::cout << "Computation time: " << (float)(duration_1 / 1000.0) << " sec" << std::endl;

    // Generamos el nombre del fichero
    stringstream ss;
    ss << argv[2];
    string str = ss.str();
    string ficheroGuardar = str;

    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Execution time: " << (float)(duration / 1000.0) << " sec" << std::endl;

    rc = MPI_Finalize();
}