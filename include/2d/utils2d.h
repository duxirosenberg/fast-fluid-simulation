#ifndef UTILS_2D_H
#define UTILS_2D_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lbm2d_struct.h"

// Shared function pointer type for LBM simulation functions
typedef void(*comp_func)(double *, int *, int, int, int);
//function pointers for timing
typedef void(*f_structs)(struct LBM*);


#ifdef __cplusplus
extern "C" {
#endif
// Function prototypes for registering and performing simulations
void add_function(comp_func f, const char* name, int flops);

#ifdef __cplusplus
}
#endif


// Constants

#define CYCLES_REQUIRED 1e8
#define REP 5
#define EPS 1e-3
#define Nt 1000
#define NL 9
#define rho0 100
#define taus 0.6
#define tol 1e-6
#define tau 0.6

static const int cxs[NL] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
static const int cys[NL] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
static const double weights[NL] = {4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36};


// function to verify two results are the same // up to 6 decimal places
int equal(double *arr1, double *arr2, int size) {
    for (int i = 0; i < size; i++) {
        double diff = arr1[i] - arr2[i];
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


void initialize_F(double *F, int Nx, int Ny){
    double* rho = (double*)malloc(Nx * Ny * sizeof(double));

    for (int i = 0; i < Ny * Nx * NL; i++) {
        F[i] = 1.0;
    }
    // Seed and generate noise
    srand(42); // np.random.seed(42)
    for (int i = 0; i < Ny * Nx * NL; i++) {
        F[i] += 0.01 * (((double)rand() / RAND_MAX) * 2.0 - 1.0); // Generate numbers between -1 and 1
    }
    // Apply the initial conditions
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            int index = y * Nx * NL + x * NL + 3;
            F[index] += 2 * (1 + 0.2 * cos(2 * M_PI * x / Nx * 4));
        }
    }
    // Compute rho
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            int rho_index = y * Nx + x;
            rho[rho_index] = 0; // Initialize rho
            for (int l = 0; l < NL; l++) {
                int index = y * Nx * NL + x * NL + l;
                rho[rho_index] += F[index];
            }
        }
    }
    // Apply the adjustment for each layer based on rho
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            for (int l = 0; l < NL; l++) {
                int index = y * Nx * NL + x * NL + l;
                int rho_index = y * Nx + x;
                F[index] *= rho0 / rho[rho_index];
            }
        }
    }
}    


void initialize_cylinder(int *cylinder, int Nx, int Ny){
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            if ((j-Nx/4)*(j-Nx/4) + (i-Ny/2)*(i-Ny/2) < (Ny/4)*(Ny/4)) {
                cylinder[j + i * Nx] = 1; // Inside the cylinder
            } else {
                cylinder[j + i * Nx] = 0; // Outside the cylinder
            }
        }
    }
}



#endif // UTILS_2D_H


