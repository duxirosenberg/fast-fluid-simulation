#include "2D_collision.h"
#include "lbm2d_struct.h"
            
#include "utils2d.h"
#include <stdio.h>
#include <stdlib.h>


/*---------------------------------BASELINE  -Do not modify!!!!!!! -------------------------------------------*/

void collision(struct LBM* S) {
    /*
    Function that performs the collision step of the LBM simulation.
    */

    //calculate fluid variables
    for (int y = 0; y < S->Ny; ++y) {
            for (int x = 0; x < S->Nx; ++x) {
                double sum_f = 0;
                double sum_fx = 0;
                double sum_fy = 0;
                for (int l = 0; l < NL; ++l) {
                    int index = (y * S->Nx + x) * NL + l;
                    double f_val = S->F[index];  
                    sum_f += f_val;
                    sum_fx += f_val * cxs[l];
                    sum_fy += f_val * cys[l];
                }
                int cell_index = y * S->Nx + x;
                S->rho[cell_index] = sum_f;
                S->ux[cell_index] = sum_fx / sum_f;
                S->uy[cell_index] = sum_fy / sum_f;
            }
    }
    //collision
    for (int y = 0; y < S->Ny; ++y) {
        for (int x = 0; x < S->Nx; ++x) {
            int idx = y * S->Nx + x; // Index for the current cell
            for (int l = 0; l < NL; ++l) {
                double cx = cxs[l];
                double cy = cys[l];
                double w = weights[l];
                double cu = 3 * (cx * S->ux[idx] + cy * S->uy[idx]);
                double uu = S->ux[idx] * S->ux[idx] + S->uy[idx] * S->uy[idx];
                S->Feq[idx * NL + l] = S->rho[idx] * w * (1 + cu + 0.5 * (cu * cu) - 1.5 * uu);
            }
        }
    }
    for (int i = 0; i < S->Nx * S->Ny * NL; ++i) {
        S->F[i] += -(1.0 / tau) * (S->F[i] - S->Feq[i]);
    }
}



/*--------------------------------- Add optimized drift functions here -------------------------------------------*/

void opt01(struct LBM* S) {
    collision(S);
}




void register_collision_functions(){
    add_collision_func(collision, "Baseline Collision",0);
    // Add optimized drift functions BELOW:

    //ONLY MODIFY HERE:
    add_collision_func(opt01, "Optimized Collision 01",0);
}