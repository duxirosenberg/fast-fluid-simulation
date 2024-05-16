
#ifndef CMDLINE_LBM_MOMENTUM_H
#define CMDLINE_LBM_MOMENTUM_H
#include "timing.h"

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
// Doubles: // particle_dist: nXYZ*q;
            // velocity_field: nXYZ*3  write
            // denssity field nXYZ,    write  
// Ints    // dirs: 3xq
void momentum_baseline(struct LBMarrays* S);

void momentumO1(struct LBMarrays* S);
void momentumO2(struct LBMarrays* S);

void momentum_arrays(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions);


static struct ops momentum_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (3 + 4 * S->direction_size),
            val * (10 + 14 * S->direction_size),
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}

static struct ops momentum_O1_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (4 * S->direction_size),
            3 + val * (5 + S->direction_size) + S->nX * S->nY,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}

static struct ops momentum_O2_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (4 * S->direction_size),
            3 + val * (3 + 4*S->direction_size) + S->nX * S->nY,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}


static void register_momentum_functions() {
    add_momentum_array_func(&momentum_arrays, &momentum_baseline_flops, "Momentum - Arrays Bl");
    add_momentum_struct_func(&momentum_baseline, &momentum_baseline_flops, "Momentum - Structs Bl");
    add_momentum_struct_func(&momentumO1, &momentum_O1_flops, "Momentum - 01");
    add_momentum_struct_func(&momentumO2, &momentum_O2_flops, "Momentum - only some optimizations (should be worse than O1) ");
}

#endif //CMDLINE_LBM_MOMENTUM_H