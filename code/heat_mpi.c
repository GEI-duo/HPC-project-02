#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BMP_HEADER_SIZE 54
#define ALPHA 0.01 // Thermal diffusivity
#define L 0.2      // Length (m) of the square domain
#define DX 0.02    // grid spacing in x-direction
#define DY 0.02    // grid spacing in y-direction
#define DT 0.0005  // Time step
#define T 1500     // Temperature on ºk of the heat source

#include <omp.h>
#include <mpi.h>

void initialize_grid(double *grid, int nx, int ny, int rank, int size, int local_nx)
{
    int base_rows = nx / size;
    int extra = nx % size;
    int row_offset = rank * base_rows + (rank < extra ? rank : extra);
    int i, j, global_i;
    double value;

#pragma omp parallel for collapse(2) private(i, j, global_i, value)
    for (i = 1; i <= local_nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            global_i = row_offset + i - 1;
            value = 0.0;

            if (global_i != 0 && global_i != nx - 1 && j != 0 && j != ny - 1)
            {
                if (global_i == j || global_i == nx - 1 - j)
                {
                    value = T;
                }
            }

            grid[i * ny + j] = value;
        }
    }
}


void solve_heat_equation(double *grid, double *new_grid, int steps, double r, int ny, int local_nx, int rank, int size)
{
    double *temp;
    int i, j, step, global_i;
    int is_first = (rank == 0);
    int is_last = (rank == size - 1);
    int idx;

    for (step = 1; step < steps; step++)
    {
        // Exchange ghost rows
        if (!is_first)
        {
            MPI_Sendrecv(&grid[1 * ny], ny, MPI_DOUBLE, rank - 1, 0,
                         &grid[0 * ny], ny, MPI_DOUBLE, rank - 1, 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (!is_last)
        {
            MPI_Sendrecv(&grid[local_nx * ny], ny, MPI_DOUBLE, rank + 1, 1,
                         &grid[(local_nx + 1) * ny], ny, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

#pragma omp parallel for private(i, j) collapse(2)
        for (i = 1; i <= local_nx; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                idx = i * ny + j;
                new_grid[idx] = grid[idx]
                                + r * (grid[idx + ny] + grid[idx - ny] - 2 * grid[idx])
                                + r * (grid[idx + 1] + grid[idx - 1] - 2 * grid[idx]);
            }
        }

        // Set Dirichlet boundaries on top and bottom rows
        if (is_first)
        {
            #pragma omp parallel for
            for (j = 0; j < ny; j++)
            {
                new_grid[1 * ny + j] = 0.0;  // First interior row, boundary above
            }
        }
        if (is_last)
        {
            #pragma omp parallel for
            for (j = 0; j < ny; j++)
            {
                new_grid[local_nx * ny + j] = 0.0;  // Last interior row, boundary below
            }
        }

        temp = grid;
        grid = new_grid;
        new_grid = temp;
    }
}

void write_bmp_header(FILE *file, int width, int height)
{
    unsigned char header[BMP_HEADER_SIZE] = {0};

    int file_size = BMP_HEADER_SIZE + 3 * width * height;
    header[0] = 'B';
    header[1] = 'M';
    header[2] = file_size & 0xFF;
    header[3] = (file_size >> 8) & 0xFF;
    header[4] = (file_size >> 16) & 0xFF;
    header[5] = (file_size >> 24) & 0xFF;
    header[10] = BMP_HEADER_SIZE;

    header[14] = 40; // Info header size
    header[18] = width & 0xFF;
    header[19] = (width >> 8) & 0xFF;
    header[20] = (width >> 16) & 0xFF;
    header[21] = (width >> 24) & 0xFF;
    header[22] = height & 0xFF;
    header[23] = (height >> 8) & 0xFF;
    header[24] = (height >> 16) & 0xFF;
    header[25] = (height >> 24) & 0xFF;
    header[26] = 1;  // Planes
    header[28] = 24; // Bits per pixel

    fwrite(header, 1, BMP_HEADER_SIZE, file);
}

void get_color(double value, unsigned char *r, unsigned char *g, unsigned char *b)
{

    if (value >= 500.0)
    {
        *r = 255;
        *g = 0;
        *b = 0; // Red
    }
    else if (value >= 100.0)
    {
        *r = 255;
        *g = 128;
        *b = 0; // Orange
    }
    else if (value >= 50.0)
    {
        *r = 171;
        *g = 71;
        *b = 188; // Lilac
    }
    else if (value >= 25)
    {
        *r = 255;
        *g = 255;
        *b = 0; // Yellow
    }
    else if (value >= 1)
    {
        *r = 0;
        *g = 0;
        *b = 255; // Blue
    }
    else if (value >= 0.1)
    {
        *r = 5;
        *g = 248;
        *b = 252; // Cyan
    }
    else
    {
        *r = 255;
        *g = 255;
        *b = 255; // white
    }
}

void write_grid(FILE *file, double *grid, int nx, int ny)
{
    int row_stride = ny * 3;
    int padding = (4 - (row_stride % 4)) % 4;
    int padded_row_size = row_stride + padding;
    int total_size = nx * padded_row_size;

    unsigned char *pixel_data = malloc(total_size);
    if (!pixel_data)
    {
        fprintf(stderr, "Failed to allocate memory\n");
        exit(1);
    }

    int i, j, p;

#pragma omp parallel for private(i, j, p)
    for (i = 0; i < nx; i++)
    {
        int row_index = nx - 1 - i;
        int row_offset = i * padded_row_size;

        for (j = 0; j < ny; j++)
        {
            unsigned char r, g, b;
            get_color(grid[row_index * ny + j], &r, &g, &b);
            int pixel_offset = row_offset + j * 3;
            pixel_data[pixel_offset + 0] = b;
            pixel_data[pixel_offset + 1] = g;
            pixel_data[pixel_offset + 2] = r;
        }

        for (p = 0; p < padding; p++)
        {
            pixel_data[row_offset + row_stride + p] = 0;
        }
    }

    fwrite(pixel_data, 1, total_size, file);
    free(pixel_data);
}

// Main function
int main(int argc, char *argv[])
{
    int nx, ny, steps;
    double start_time, end_time;
    char car;
    int rank, size;

    // We will split into a GRID and assign each worker a cell of the grid
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4)
    {
        printf("Usage: heat_mpi <size> <steps> <name_output_file>\n");
        MPI_Finalize();
        return 1;
    }

    double r = ALPHA * DT / (DX * DY);
    nx = ny = atoi(argv[1]);
    steps = atoi(argv[2]);

    if (rank == 0)
        start_time = MPI_Wtime();

    // Local size of the process
    int base_rows = nx / size;
    int remainder = nx % size;
    int local_nx = base_rows + (rank < remainder ? 1 : 0);

    // Allocate grid with 2 extra rows for ghost rows
    double *grid = (double *)calloc((local_nx + 2) * ny, sizeof(double));
    double *new_grid = (double *)calloc((local_nx + 2) * ny, sizeof(double));

    initialize_grid(grid, nx, ny, rank, size, local_nx);
    solve_heat_equation(grid, new_grid, steps, r, ny, local_nx, rank, size);

    double *full_grid = NULL;
    int *recvcounts = NULL, *displs = NULL;

    if (rank == 0)
    {
        full_grid = malloc(nx * ny * sizeof(double));
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));

        int offset = 0;
        int i;
        for (i = 0; i < size; i++)
        {
            int rows = base_rows + (i < remainder ? 1 : 0);
            recvcounts[i] = rows * ny;
            displs[i] = offset;
            offset += rows * ny;
        }
    }

    MPI_Gatherv(&grid[1 * ny], local_nx * ny, MPI_DOUBLE,
                full_grid, recvcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        FILE *file = fopen(argv[3], "wb");
        if (!file)
        {
            fprintf(stderr, "Failed to open output file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        write_bmp_header(file, nx, ny);
        write_grid(file, full_grid, nx, ny);
        fclose(file);

        end_time = MPI_Wtime();
        printf("Execution Time = %fs for %dx%d grid and %d steps\n", end_time - start_time, nx, ny, steps);

        free(full_grid);
        free(recvcounts);
        free(displs);
    }

    // Free allocated memory
    free(grid);
    free(new_grid);
    MPI_Finalize();

    return 0;
}
