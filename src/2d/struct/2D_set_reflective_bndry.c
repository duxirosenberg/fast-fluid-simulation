#include "lbm2d_struct.h"
#include "utils2d.h"
#include <stdio.h>
#include <stdlib.h>


/*---------------------------------BASELINE  -Do not modify!!!!!!! -------------------------------------------*/

void set_reflective_boundaries(struct LBM* S) {
    /*
    Function that sets the reflective boundaries of the LBM simulation.
    */
    int rearrange_idx[9] = {0, 5, 6, 7, 8, 1, 2, 3, 4}; // The rearrangement pattern
    int idx = 0;
    for (int y = 0; y < S->Ny; ++y) {
        for (int x = 0; x < S->Nx; ++x) {
            if (S->cylinder[y * S->Nx + x]) { // If cylinder at (y, x) is non-zero
                int base_idx = (y * S->Nx + x) * 9; // Base index for the (y, x) slice in F_flat
                double rearranged[9];
                for (int i = 0; i < 9; ++i) {
                    rearranged[i] = S->F[base_idx + rearrange_idx[i]];
                }
                for (int i = 0; i < 9; ++i) {
                    S->bndryF[idx + i] = rearranged[i];
                }
                idx += 9;
            }
        }
    }
}


/*---------------------------------Add optimized drift functions here -------------------------------------------*/

void opt01(struct LBM* S) {
    set_reflective_boundaries(S);
}






void register_set_reflective_boundaries_functions(){
    add_reflective_boundaries_func(set_reflective_boundaries, "Baseline Set Reflective Boundaries",0);
    // Add optimized drift functions BELOW:

    //ONLY MODIFY HERE:
    add_reflective_boundaries_func(opt01, "Optimized Set Reflective Boundaries 01",0);
}