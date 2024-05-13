#include "utils2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void drift(double *F, double *temp_F, int Nx, int Ny){
    /*
    Function that performs the drift step of the LBM simulation.
    */

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
}


void set_reflective_boundaries(double *F, int *cylinder, double *bndryF, int Nx, int Ny){
    /*
    Function that sets the reflective boundaries of the LBM simulation.
    */
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
}


void collision(double *F, double *Feq, double *rho, double *ux, double *uy, int Nx, int Ny){
    /*
    Function that performs the collision step of the LBM simulation.
    */
    //calculate fluid variables
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
    //collision
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
    for (int i = 0; i < Nx * Ny * NL; ++i) {
        F[i] += -(1.0 / tau) * (F[i] - Feq[i]);
    }

}

void apply_boundary(double *F, int *cylinder, double *bndryF, int Nx, int Ny){
    /*
    Function that applies the boundary conditions of the LBM simulation.
    */
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
}


void lbm_baseline(double *F, int *cylinder, int Nx, int Ny, int n_steps){
    /*
    Function that takes an initial particle distribution F and performs n_steps updates of the LBM simulation
    */

    // Allocate memory
    double* rho = (double*)calloc(Nx * Ny,sizeof(double));
    double* ux = (double*)calloc(Nx * Ny , sizeof(double));
    double* uy = (double*)calloc(Nx * Ny , sizeof(double));
    double* Feq = (double*)calloc(Nx * Ny * NL , sizeof(double));
    //TODO: bndryF does not need Nx * Ny * 9, it needs only the number of cylinder cells * 9
    double* bndryF = (double*)calloc(Nx * Ny * NL, sizeof(double)); // Allocate memory for bndryF
    double* temp_F = (double*)malloc(Nx * Ny * NL* sizeof(double)); 



    for (int step = 0; step < n_steps; step++) {
                   
        //DRIFT:
        drift(F, temp_F, Nx, Ny);

        // SET REFLECTIVE BOUNDARIES
        set_reflective_boundaries(F, cylinder, bndryF, Nx, Ny);

        //3. CALCULATE FLUID VARIABLES and
        // APPLY COLLISION
        collision(F, Feq, rho, ux, uy, Nx, Ny);

        // APPLY BOUNDARY
        apply_boundary(F, cylinder, bndryF, Nx, Ny);
    }

    // Free memory
    free(temp_F);
    free(bndryF);
    free(rho);
    free(ux);
    free(uy);
    free(Feq);

}


void lbm_opt01(double *F, int *cylinder, int Nx, int Ny, int n_steps){
    lbm_baseline(F, cylinder, Nx, Ny, n_steps);
}


// Function to register LBM functions
void register_functions(){
    //Baseline Runs Automatically 

    // add more functions here using
    // add_function(function_name, "Description", #flops per Nx*Ny);
    add_function(lbm_opt01, "First optimization", 808);
}


