
#ifndef CMDLINE_LBM_COUETTE_H
#define CMDLINE_LBM_COUETTE_H
#include "timing.h"

// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
void stream_couette_baseline(struct LBMarrays* S);

void stream_couette_arrays(int nX, int nY, int nZ, int direction_size, double c_s,
                           double* previous_particle_distributions,
                           double* particle_distributions,
                           const int* directions,
                           const double* weights,
                           const int* reverse_indexes);

static struct ops stream_couette_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    struct ops ops = {
            5 * val / 3,
            val * (39 * S->nY + 32 / 3)
    };
    return ops;
}

static void register_stream_couette_functions() {
    add_stream_couette_array_func(&stream_couette_arrays, &stream_couette_baseline_flops, "Stream Couette Baseline - Arrays");
    add_stream_couette_struct_func(&stream_couette_baseline, &stream_couette_baseline_flops, "Stream Couette Baseline - Structs");
}

#endif //CMDLINE_LBM_COUETTE_H
