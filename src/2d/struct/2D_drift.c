#include "lbm2d_struct.h"
#include "utils2d.h"
#include <stdio.h>
#include <stdlib.h>


/*---------------------------------BASELINE-------------------------------------------*/
void drift(struct LBM* S) {
 /*
    Function that performs the drift step of the LBM simulation.
    */
    // Drift in the x direction 
    for (int l = 0; l < NL; l++) {
        for (int y = 0; y < S->Ny; y++) {
            for (int x = 0; x < S->Nx; x++) {
                int source_index = y * S->Nx * NL + x * NL + l;
                int shifted_x = ((x + cxs[l]) % S->Nx + S->Nx) % S->Nx; // The +Nx % Nx is to handle negative numbers
                int target_index = y * S->Nx * NL + shifted_x * NL + l;
                S->temp_F[target_index] = S->F[source_index];
            }
        }
    }
    for (int i = 0; i < S->Nx * S->Ny * NL; ++i) {
        S->F[i] = S->temp_F[i];
    }
    // Drift in the y direction 
    for (int l = 0; l < NL; l++) {
        for (int x = 0; x < S->Nx; x++) {
            for (int y = 0; y < S->Ny; y++) {
                int source_index = y * S->Nx * NL + x * NL + l;
                int shifted_y = ((y + cys[l]) % S->Ny + S->Ny) % S->Ny; // The +Ny % Ny is to handle negative numbers
                int target_index = shifted_y * S->Nx * NL + x * NL + l;
                S->temp_F[target_index] = S->F[source_index];
            }
        }
    }
    for (int i = 0; i < S->Nx * S->Ny * NL; i++) {
       S->F[i] = S->temp_F[i];
    }  
}

/*---------------------------------Add optimized drift functions here -------------------------------------------*/

void drift_opt01(struct LBM* S) {
    drift(S);
}


void register_apply_boundary_functions(){
    add_drift_func(drift, "Baseline Drift",0);
    // Add optimized drift functions BELOW:
    add_drift_func(drift_opt01, "Optimized Drift 01",0);
}