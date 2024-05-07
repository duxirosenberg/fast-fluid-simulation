
#ifndef CMDLINE_LBM_LEES_EDWARDS_H
#define CMDLINE_LBM_LEES_EDWARDS_H
#include "timing.h"

// Flops: 4 + (nX * nZ * q / 3) * 233
// Intops: 1 + nX * nZ * q * (32 * nY + 106 / 3)
void stream_lees_edwards_baseline(struct LBMarrays* S, int time);

void stream_lees_edwards_arrays(int nX, int nY, int nZ, int direction_size, int time, double gamma_dot, double c_s,
                                double* density_field,
                                double* velocity_field,
                                double* previous_particle_distributions,
                                double* particle_distributions,
                                const int* directions,
                                const double* weights
);

static struct ops stream_lees_edwards_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    struct ops ops = {
            4 + (val / 3 * 233),
            1 + (val * (32 * S->nY + 106 / 3))
    };
    return ops;
}

static void register_stream_lees_edwards_functions() {
    add_stream_lees_edwards_array_func(&stream_lees_edwards_arrays, &stream_lees_edwards_baseline_flops, "Stream Lees Edwards Baseline - Arrays");
    add_stream_lees_edwards_struct_func(&stream_lees_edwards_baseline, &stream_lees_edwards_baseline_flops, "Stream Lees Edwards Baseline - Structs");
}

#endif //CMDLINE_LBM_LEES_EDWARDS_H