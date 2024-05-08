
#ifndef CMDLINE_LBM_MOMENTUM_H
#define CMDLINE_LBM_MOMENTUM_H
#include "timing.h"

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
// Doubles: // particle_dist: nXYZ*q;
            // velocity_field: nXYZ*3
            // denssity field nXYZ, 
// Ints    // dirs: 3xq
void momentum_baseline(struct LBMarrays* S);

void momentum_arrays(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions);

static struct ops momentum_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (3 + 4 * S->direction_size),
            val * (10 + 14 * S->direction_size),
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int))
    };
    return ops;
}

static void register_momentum_functions() {
    add_momentum_array_func(&momentum_arrays, &momentum_baseline_flops, "Momentum Baseline - Arrays");
    add_momentum_struct_func(&momentum_baseline, &momentum_baseline_flops, "Momentum Baseline - Structs");
}

#endif //CMDLINE_LBM_MOMENTUM_H
