#include "timing.h"

#ifndef CMDLINE_LBM_COUETTE_H
#define CMDLINE_LBM_COUETTE_H

void stream_couette_baseline(struct LBMarrays* S);

void stream_couette_arrays(int nX, int nY, int nZ, int direction_size, double c_s,
                           double* previous_particle_distributions,
                           double* particle_distributions,
                           const int* directions,
                           const double* weights,
                           const int* reverse_indexes);

static void register_stream_couette_functions() {
    // TODO johannes: flops
    add_stream_couette_array_func(&stream_couette_arrays, "Stream Couette Baseline - Arrays", 10);
    add_stream_couette_struct_func(&stream_couette_baseline, "Stream Couette Baseline - Structs", 10);
}

#endif //CMDLINE_LBM_COUETTE_H
