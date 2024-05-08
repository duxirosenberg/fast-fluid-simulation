
#ifndef CMDLINE_LBM_PERIODIC_H
#define CMDLINE_LBM_PERIODIC_H
#include "timing.h"

// Flops: 0
// Intops: xN * nY * nZ * q * 32
// Doubles: // particle_dist: nXYZ*q;
            // pev_particle_dist: nXYZ*q;
// Ints     // dirs: 3xq
void stream_periodic_baseline(struct LBMarrays* S, int time);

void stream_periodic_arrays(int nX, int nY, int nZ, int direction_size,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions
);

static struct ops stream_periodic_baseline_flops(struct LBMarrays* S) {

    long val = S->nX * S->nY * S->nZ* S->direction_size;
    struct ops ops = {
            0,
            val * 32,
            val*2*(int)(sizeof(double)) + 3*S->direction_size*(int)(sizeof(int))
    };
    return ops;
}

static void register_stream_periodic_functions() {
    add_stream_periodic_array_func(&stream_periodic_arrays, &stream_periodic_baseline_flops, "Stream Periodic Baseline - Arrays");
    add_stream_periodic_struct_func(&stream_periodic_baseline, &stream_periodic_baseline_flops, "Stream Periodic Baseline - Structs");
}

#endif //CMDLINE_LBM_PERIODIC_H
