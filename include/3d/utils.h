#ifndef COMMON_H
#define COMMON_H
#include <cmath>
#include <cstdio>
#include <cstdlib>

extern "C" {
    #include "LBM.h"
    #include "LBMFunctions.h"
    #include "momentum.h"
    #include "collision.h"
    #include "couette.h"
    #include "periodic.h"
    #include "lees_edwards.h"
    #include "timing.h"
}

// Constants
#define tol 1e-10


static const int D2Q9_DIRECTIONS[] = { 1,0,0, 0,1,0, -1,0,0, 0,-1,0, 1,1,0, -1,1,0, -1,-1,0, 1,-1,0, 0,0,0};
static const double D2Q9_WEIGHTS[] = {1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,4.0/9.0};

static const int D3Q15_DIRECTIONS[] = { 0,0,0, 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1, 1,1,1, -1,-1,-1, 1,1,-1, -1,-1,1, 1,-1,1, -1,1,-1, -1,1,1, 1,-1,-1 };
static const double D3Q15_WEIGHTS[] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0  };

static const int D3Q27_DIRECTIONS[] = { 0,0,0,  1,0,0, -1,0,0,0,1,0, 0,-1,0, 0,0,1,
                                 0,0,-1, 1,1,0, -1,-1,0, 1,0,1, -1,0,-1, 0,1,1,
                                 0,-1,-1, 1,-1,0, -1,1,0, 1,0,-1, -1,0,1, 0,1,-1,
                                 0,-1,1, 1,1,1, -1,-1,-1, 1,1,-1, -1,-1,1, 1,-1,1,
                                 -1,1,-1, -1,1,1, 1,-1,-1};
static const double D3Q27_WEIGHTS[] = {8.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0, 1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0,1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};


// function to verify two results are the same // up to 6 decimal places
int equal(double *arr1, double *arr2, int size) {
    for (int i = 0; i < size; i++) {
        double diff = fabs(arr1[i] - arr2[i]);
        if (diff > tol) {
            return 0;
        }
    }
    return 1;
}


void read_csv_to_array(const char* filename, double* array, int size) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Could not open file %s\n", filename);
        return;
    }
    // Read values from CSV file and store them into the array
    for (int i = 0; i < size; ++i) {
        if (fscanf(file, "%lf,", &array[i]) != 1) {
            printf("Error: Failed to read value from file %s\n", filename);
            fclose(file);
            return;
        }
    }
    fclose(file);
}


void write_array_to_csv(const char* filename, double* array, int size) {
    FILE* file = fopen(filename, "w");
    for (int i = 0; i < size;i++){
        fprintf(file, "%lf\n", array[i]);
    }
    fclose(file);
}

void initialize(int nX, int nY, int nZ, int direction_size,
                double* density_field,
                double* velocity_field,
                double* previous_particle_distributions,
                double* particle_distributions,
                const int* directions,
                const double* weights,
                int* reverse_indexes
) {


    for(int i = 0; i < direction_size; i++) {
        for(int j = 0; j < direction_size; j++) {
            if(directions[3 * i] == -directions[3 * j] && directions[3 * i + 1] == -directions[3 * j + 1] && directions[3 * i + 2] == -directions[3 * j + 2]) {
                reverse_indexes[i] = j;
            }
        }
    }

    int box_flatten_length = nX * nY * nZ;

    for(int i = 0; i < box_flatten_length; i++) {
        density_field[i] = 1;
    }

    for(int i = 0; i < 3 * box_flatten_length; i++) {
        velocity_field[i] = 0;
    }

    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            for(int z = 0; z < nZ; z++) {
                for(int i = 0; i < direction_size; i++) {
                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    previous_particle_distributions[index] = x+ weights[i];
                    particle_distributions[index] = x+weights[i];
                }
            }
        }
    }
}
#endif // COMMON_H


