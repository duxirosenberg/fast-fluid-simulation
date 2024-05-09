
#ifndef CMDLINE_LBM_COLLISION_H
#define CMDLINE_LBM_COLLISION_H
#include "timing.h"

// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
// Doubles: 
            // prev_particle_dist: nXYZ*q;    write
            // velocity_field: nXYZ*3
            // denssity field nXYZ, 
            // weights q
//  Ints:   // dirs: 3xq
void collision_arrays(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights);

void collision_baseline(struct LBMarrays* S);

static struct ops collision_baseline_flops(struct LBMarrays* S) {
    long val0 = S->nX * S->nY * S->nZ;
    long val = val0 * S->direction_size;
    struct ops ops = {
            2 + val * 29,
            val * 23,
            val0*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(3*sizeof(int) +sizeof(double)),
            val*(int)(sizeof(double))
    
    };
    return ops;
}

static void register_collision_functions() {
    add_collision_array_func(&collision_arrays, &collision_baseline_flops, "Collision - Arrays Bl");
    add_collision_struct_func(&collision_baseline, &collision_baseline_flops, "Collision - Structs Bl");
}

#endif //CMDLINE_LBM_COLLISION_H
