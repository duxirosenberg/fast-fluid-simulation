
#ifndef CMDLINE_LBM_PERIODIC_H
#define CMDLINE_LBM_PERIODIC_H
#include "timing.h"

// Flops: 0
// Intops: xN * nY * nZ * q * 32
void stream_periodic_baseline(struct LBMarrays* S, int time);

void stream_periodic_arrays(int nX, int nY, int nZ, int direction_size,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions
);

static struct ops stream_periodic_baseline_flops(struct LBMarrays* S) {
    struct ops ops = {
            0,
            S->nX * S->nY * S->nZ * S->direction_size * 32
    };
    return ops;
}

static void register_stream_periodic_functions() {
    add_stream_periodic_array_func(&stream_periodic_arrays, &stream_periodic_baseline_flops, "Stream Periodic Baseline - Arrays");
    add_stream_periodic_struct_func(&stream_periodic_baseline, &stream_periodic_baseline_flops, "Stream Periodic Baseline - Structs");
}

#endif //CMDLINE_LBM_PERIODIC_H
