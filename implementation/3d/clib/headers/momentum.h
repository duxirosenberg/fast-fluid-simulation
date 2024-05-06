#include "timing.h"

#ifndef CMDLINE_LBM_MOMENTUM_H
#define CMDLINE_LBM_MOMENTUM_H

void momentum_baseline(struct LBMarrays* S);

void momentum_arrays(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions);

static void register_momentum_functions() {
    // TODO johannes: flops
    add_momentum_array_func(&momentum_arrays, "Momentum Baseline - Arrays", 10);
    add_momentum_struct_func(&momentum_baseline, "Momentum Baseline - Structs", 10);
}

#endif //CMDLINE_LBM_MOMENTUM_H
