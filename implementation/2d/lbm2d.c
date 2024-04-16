#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define DEBUG 0

int cxs[NL] = {0, 0, 1, 1, 1, 0,-1,-1,-1};
int cys[NL] = {0, 1, 1, 0,-1,-1,-1, 0, 1};
double weights[NL] = {4.0/9,1.0/9,1.0/36,1.0/9,1.0/36,1.0/9,1.0/36,1.0/9,1.0/36}; //sums to 1

void lbm_baseline(double *F, int *cylinder, int Nx, int Ny, int n_steps){
    /*
    Function that takes an initial particle distribution F and performs n_steps updates of the LBM simulation
    */


    for (int step = 0; step < n_steps; step++) {
        /*
        if(step % 100 == 0){
            printf("C Baseline: timestep %d/%d\n", step, n_steps);
        }
        */
        
        // Allocate memory
        double* rho = (double*)calloc(Nx * Ny,sizeof(double));
        double* ux = (double*)calloc(Nx * Ny , sizeof(double));
        double* uy = (double*)calloc(Nx * Ny , sizeof(double));
        double* Feq = (double*)calloc(Nx * Ny * NL , sizeof(double));
        //TODO: bndryF does not need Nx * Ny * 9, it needs only the number of cylinder cells * 9
        double* bndryF = (double*)calloc(Nx * Ny * 9, sizeof(double)); // Allocate memory for bndryF

            
        //DRIFT:
        double* temp_F = (double*)malloc(Nx * Ny * NL* sizeof(double)); 
        // Drift in the x direction 
        for (int l = 0; l < NL; l++) {
            for (int y = 0; y < Ny; y++) {
                for (int x = 0; x < Nx; x++) {
                    int source_index = y * Nx * NL + x * NL + l;
                    int shifted_x = ((x + cxs[l]) % Nx + Nx) % Nx; // The +Nx % Nx is to handle negative numbers
                    int target_index = y * Nx * NL + shifted_x * NL + l;
                    temp_F[target_index] = F[source_index];
                }
            }
        }
        for (int i = 0; i < Nx * Ny * NL; ++i) {
            F[i] = temp_F[i];
        }
        // Drift in the y direction 
        for (int l = 0; l < NL; l++) {
            for (int x = 0; x < Nx; x++) {
                for (int y = 0; y < Ny; y++) {
                    int source_index = y * Nx * NL + x * NL + l;
                    int shifted_y = ((y + cys[l]) % Ny + Ny) % Ny; // The +Ny % Ny is to handle negative numbers
                    int target_index = shifted_y * Nx * NL + x * NL + l;
                    temp_F[target_index] = F[source_index];
                }
            }
        }
        for (int i = 0; i < Nx * Ny * NL; i++) {
            F[i] = temp_F[i];
        }
        free(temp_F);


        if (DEBUG) {
            write_array_to_csv("../py_baseline/F_c_drift.csv", F, Nx * Ny * NL);
        }


        // SET REFLECTIVE BOUNDARIES
        int rearrange_idx[9] = {0, 5, 6, 7, 8, 1, 2, 3, 4}; // The rearrangement pattern
        int idx = 0;
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                if (cylinder[y * Nx + x]) { // If cylinder at (y, x) is non-zero
                    int base_idx = (y * Nx + x) * 9; // Base index for the (y, x) slice in F_flat
                    double rearranged[9];
                    for (int i = 0; i < 9; ++i) {
                        rearranged[i] = F[base_idx + rearrange_idx[i]];
                    }
                    for (int i = 0; i < 9; ++i) {
                        bndryF[idx + i] = rearranged[i];
                    }
                    idx += 9;
                }
            }
        }

        if (DEBUG) {
            write_array_to_csv("../py_baseline/bndryF_c.csv", bndryF, Nx*Ny);
        }

        // CALCULATE FLUID VARIABLES
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                double sum_f = 0;
                double sum_fx = 0;
                double sum_fy = 0;
                for (int l = 0; l < NL; ++l) {
                    int index = (y * Nx + x) * NL + l;
                    double f_val = F[index];  
                    sum_f += f_val;
                    sum_fx += f_val * cxs[l];
                    sum_fy += f_val * cys[l];
                }
                int cell_index = y * Nx + x;
                rho[cell_index] = sum_f;
                ux[cell_index] = sum_fx / sum_f;
                uy[cell_index] = sum_fy / sum_f;
            }
        }
        if (DEBUG) {
            write_array_to_csv("../py_baseline/rho_c.csv", rho, Nx * Ny);
            write_array_to_csv("../py_baseline/ux_c.csv", ux, Nx * Ny);
            write_array_to_csv("../py_baseline/uy_c.csv", uy, Nx * Ny);
        }


        // APPLY COLLISION
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                int idx = y * Nx + x; // Index for the current cell
                for (int l = 0; l < NL; ++l) {
                    double cx = cxs[l];
                    double cy = cys[l];
                    double w = weights[l];
                    double cu = 3 * (cx * ux[idx] + cy * uy[idx]);
                    double uu = ux[idx] * ux[idx] + uy[idx] * uy[idx];
                    Feq[idx * NL + l] = rho[idx] * w * (1 + cu + 0.5 * (cu * cu) - 1.5 * uu);
                }
            }
        }

        if (DEBUG) {
            write_array_to_csv("../py_baseline/Feq_c.csv", Feq, Nx * Ny * NL);
        }


        for (int i = 0; i < Nx * Ny * NL; ++i) {
            F[i] += -(1.0 / tau) * (F[i] - Feq[i]);
        }


        // APPLY BOUNDARY
        int bndry_idx = 0;
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                if (cylinder[y * Nx + x]) {
                    for (int l = 0; l < NL; ++l) {
                        int F_index = (y * Nx * NL) + (x * NL) + l;
                        // Update F with the corresponding bndryF value
                        F[F_index] = bndryF[bndry_idx];
                        bndry_idx++; // Move to the next element in bndryF for the next iteration
                    }
                }
            }
        }

        if (DEBUG) {
            write_array_to_csv("../py_baseline/F_c.csv", F, Nx * Ny * NL);
        }

        // Free memory
        free(bndryF);
        free(rho);
        free(ux);
        free(uy);
        free(Feq);

    }

}

void lbm_opt01(double *F, int *cylinder, int Nx, int Ny, int n_steps){
    lbm_baseline(F, cylinder, Nx, Ny, n_steps);
}


// Function to register LBM functions
void register_functions(){
    //Baseline Runs Automatically 

    // add more functions here using
    // add_function(function_name, "Description", 1);
    
    add_function(lbm_opt01, "First optimization: ", 1);

}


