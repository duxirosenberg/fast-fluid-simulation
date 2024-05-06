
#ifndef CMDLINE_LBM_COLLISION_H
#define CMDLINE_LBM_COLLISION_H
#include "timing.h"

// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
void collision_arrays(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights);

void collision_baseline(struct LBMarrays* S);

static struct ops collision_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ * S->direction_size;
    struct ops ops = {
            2 + val * 29,
            val * 23
    };
    return ops;
}

static void register_collision_functions() {
    add_collision_array_func(&collision_arrays, &collision_baseline_flops, "Collision Baseline - Arrays");
    add_collision_struct_func(&collision_baseline, &collision_baseline_flops, "Collision Baseline - Structs");
}

#endif //CMDLINE_LBM_COLLISION_H
