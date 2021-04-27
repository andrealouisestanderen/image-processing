#include <assert.h>
#include <cmath>
#include "png++-0.2.7/png.hpp"
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

Matrix getGaussian(int height, int width, double sigma)
{
    Matrix kernel(height, Array(width));
    double sum = 0.0;
    int i, j;

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            kernel[i][j] = exp(-(i * i + j * j) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
            sum += kernel[i][j];
        }
    }

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            kernel[i][j] /= sum;
        }
    }

    return kernel;
}

Image loadImage(const char *filename)
{
    png::image<png::rgb_pixel> image(filename);
    Image imageMatrix(3, Matrix(image.get_height(), Array(image.get_width())));

    int h, w;
    for (h = 0; h < image.get_height(); h++)
    {
        for (w = 0; w < image.get_width(); w++)
        {
            imageMatrix[0][h][w] = image[h][w].red;
            imageMatrix[1][h][w] = image[h][w].green;
            imageMatrix[2][h][w] = image[h][w].blue;
        }
    }

    return imageMatrix;
}

void saveImage(Image &image, string filename)
{
    assert(image.size() == 3);

    int height = image[0].size();
    int width = image[0][0].size();
    int x, y;

    png::image<png::rgb_pixel> imageFile(width, height);

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            imageFile[y][x].red = image[0][y][x];
            imageFile[y][x].green = image[1][y][x];
            imageFile[y][x].blue = image[2][y][x];
        }
    }
    imageFile.write(filename);
}

Image applyFilter(Image &image, Matrix &filter)
{
    assert(image.size() == 3 && filter.size() != 0);

    int height = image[0].size();
    int width = image[0][0].size();
    int filterHeight = filter.size();
    int filterWidth = filter[0].size();
    int newImageHeight = height - filterHeight + 1;
    int newImageWidth = width - filterWidth + 1;
    int d, i, j, h, w;

    Image newImage(3, Matrix(newImageHeight, Array(newImageWidth)));

    for (i = 0; i < newImageHeight; i++)
    {
        for (j = 0; j < newImageWidth; j++)
        {
            for (d = 0; d < 3; d++)
            {
                for (h = i; h < i + filterHeight; h++)
                {
                    for (w = j; w < j + filterWidth; w++)
                    {
                        newImage[d][i][j] += filter[h - i][w - j] * image[d][h][w];
                    }
                }
            }
        }
    }

    return newImage;
}

Image applyFilter(Image &image, Matrix &filter, int times)
{
    Image newImage = image;
    for (int i = 0; i < times; i++)
    {
        newImage = applyFilter(newImage, filter);
    }
    return newImage;
}

int main(int argc, char *argv[])
{

    auto t1 = std::chrono::high_resolution_clock::now();

    cout << "Loading image..." << endl;
    Image image = loadImage(argv[1]);
    cout << "Applying filter..." << endl;

    auto t1_1 = std::chrono::high_resolution_clock::now();

    //MPI PARALLELIZATION STARTS HERE
    int rank, size, tag, rc, i, j, row_start, row_stop;
    double *send, *recv, *send_back, *recv_back, data;
    MPI_Status status;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Matrix filter = getGaussian(10, 10, 50.0);

    int workers = size - 1;
    int height = image[0].size();
    int width = image[0][0].size();
    int newImageHeight = height - filter.size() + 1;
    int newImageWidth = width - filter[0].size() + 1;
    int last_rows = height % workers;
    int batch_size = (height - last_rows) / workers;

    Image newImage(3, Matrix(newImageHeight, Array(newImageWidth)));

    send = (double *)malloc(sizeof(double) * width * 3);
    recv = (double *)malloc(sizeof(double) * width * 3);
    send_back = (double *)malloc(sizeof(double) * newImageWidth * 3);
    recv_back = (double *)malloc(sizeof(double) * newImageWidth * 3);

    if (rank == 0)
    {
        row_start = 0;
        for (i = 1; i < size; i++) // EACH PROCESS DOES THE BELOW
        {
            row_stop = i * batch_size + filter.size(); // +1 ??
            if (i == size - 1)
            {
                row_stop = i * batch_size + last_rows; // +1 ??
            }
            for (int row_id = row_start; row_id < row_stop; row_id++)
            {
                for (int col_id = 0; col_id < width; col_id++)
                {
                    for (int dim_id = 0; dim_id < 3; dim_id++)
                    {
                        double data = image[dim_id][row_id][col_id];
                        send[(col_id * 3) + dim_id] = data;
                        // send = [col0d0, col0d1, col0d2, col1d0, col1d1, col1d2,...,colnd0, colnd1, colnd2]
                    }
                }
                rc = MPI_Send(send, width * 3, MPI_DOUBLE, i, row_id, MPI_COMM_WORLD);
            }
            row_start += batch_size;
        }

        row_start = 0;
        for (i = 1; i < size; i++)
        {
            row_stop = i * batch_size + 1;
            if (i == size - 1)
            {
                row_stop = i * batch_size + last_rows - filter.size() + 1;
            }
            for (int row_id = row_start; row_id < row_stop; row_id++)
            {
                rc = MPI_Recv(recv_back, newImageWidth * 3, MPI_DOUBLE, i, row_id, MPI_COMM_WORLD, &status);
                //printf("\n RECeIVED ROW: %d, FROM PROCESSOR: %d\n", row_id, i);
                for (int col_id = 0; col_id < newImageWidth * 3; col_id += 3)
                {
                    newImage[0][row_id][col_id / 3] = recv_back[col_id];
                    newImage[1][row_id][col_id / 3] = recv_back[col_id + 1];
                    newImage[2][row_id][col_id / 3] = recv_back[col_id + 2];
                }
            }
            row_start += batch_size;
        }
    }

    else
    {
        row_start = (rank - 1) * batch_size; // rank 1:0, rank 2:24, rank 3:48, rank 4:72
        row_stop = row_start + batch_size + filter.size();
        height = batch_size + filter.size();
        if (rank == size - 1)
        {
            height = batch_size + last_rows;
            row_stop = row_start + batch_size + last_rows;
        }
        Image subImage(3, Matrix(height, Array(width)));
        for (int row_id = row_start; row_id < row_stop; row_id++)
        {

            rc = MPI_Recv(recv, width * 3, MPI_DOUBLE, 0, row_id, MPI_COMM_WORLD, &status);
            //printf("\nPROC: %d, RECEIVED TAG: %d\n", rank, row_id);

            for (int col_id = 0; col_id < width * 3; col_id += 3)
            {
                //printf("\nPROC: %d, TRYING TO PLACE %d\n", rank, recv[col_id]);
                //printf("\nPROC: %d, AT ROW INDEX: %d\n", rank, row_id - (rank - 1)*batch_size);
                subImage[0][row_id - row_start][col_id / 3] = recv[col_id];
                subImage[1][row_id - row_start][col_id / 3] = recv[col_id + 1];
                subImage[2][row_id - row_start][col_id / 3] = recv[col_id + 2];
            }
        }

        Image newSubImage = applyFilter(subImage, filter);

        stringstream ss;
        ss << argv[2];
        string str = ss.str();
        string ficheroGuardar = str;

        saveImage(newSubImage, ficheroGuardar);

        int newSubHeight = newSubImage[0].size();
        int newSubWidth = newSubImage[0][0].size();
        for (int row_id = 0; row_id < newSubHeight; row_id++)
        {
            tag = (rank - 1) * batch_size + row_id;
            //printf("\n PROC: %d, SENDING BACK ORIGINAL ROW: %d, FROM INTERNAL ROW: %d\n", rank, tag, row_id);
            for (int col_id = 0; col_id < newSubWidth; col_id++)
            {
                for (int dim_id = 0; dim_id < 3; dim_id++)
                {
                    double data = newSubImage[dim_id][row_id][col_id];
                    send_back[col_id + (newSubWidth * dim_id)] = data;
                }
            }
            rc = MPI_Send(send_back, newSubWidth * 3, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
    }

    //Every slave Save the subimage as image received to check

    auto t2_1 = std::chrono::high_resolution_clock::now();
    auto duration_1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2_1 - t1_1).count();
    std::cout << "Computation time: " << (float)(duration_1 / 1000.0) << " sec" << std::endl;

    cout << "Saving image..." << endl;

    // Generamos el nombre del fichero
    stringstream ss;
    ss << argv[2];
    string str = ss.str();
    string ficheroGuardar = str;

    saveImage(newImage, ficheroGuardar);
    cout << "Done!" << endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Execution time: " << (float)(duration / 1000.0) << " sec" << std::endl;

    rc = MPI_Finalize();
}