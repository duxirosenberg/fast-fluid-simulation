#include "timing.h"

#ifndef CMDLINE_LBM_COLLISION_H
#define CMDLINE_LBM_COLLISION_H

void collision_arrays(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights);

void collision_baseline(struct LBMarrays* S);

static void register_collision_functions() {
    // TODO johannes: flops
    add_collision_array_func(&collision_arrays, "Collision Baseline - Arrays", 10);
    add_collision_struct_func(&collision_baseline, "Collision Baseline - Structs", 10);
}

#endif //CMDLINE_LBM_COLLISION_H
