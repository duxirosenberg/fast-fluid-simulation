#include "lbm2d_struct.h"
#include "utils2d.h"


#include <stdio.h>
#include <stdlib.h>

#define NL 9
/*---------------------------------BASELINE  -Do not modify!!!!!!! -------------------------------------------*/
void apply_boundary(struct LBM* S) {
    /*
    Function that applies the boundary conditions of the LBM simulation.
    */
    int bndry_idx = 0;
    for (int y = 0; y < S->Ny; ++y) {
        for (int x = 0; x < S->Nx; ++x) {
            if (S->cylinder[y * S->Nx + x]) {
                for (int l = 0; l < NL; ++l) {
                    int F_index = (y * S->Nx * NL) + (x * NL) + l;
                    // Update F with the corresponding bndryF value
                    S->F[F_index] = S->bndryF[bndry_idx];
                    bndry_idx++; // Move to the next element in bndryF for the next iteration
                }
            }
        }
    }
}


/*--------------------------------- Add optimized drift functions here -------------------------------------------*/

void opt01(struct LBM* S) {
    apply_boundary(S);
}




void register_drift_functions(){
    add_apply_boundary_func(apply_boundary, "Baseline Apply Bndry",0);
    // Add optimized drift functions BELOW:

    //ONLY MODIFY HERE:
    add_apply_boundary_func(opt01, "Optimized Apply Bndry 01",0);
}