#include "timing.h"

#ifndef CMDLINE_LBM_LEES_EDWARDS_H
#define CMDLINE_LBM_LEES_EDWARDS_H

void stream_lees_edwards_baseline(struct LBMarrays* S, int time);

void stream_lees_edwards_arrays(int nX, int nY, int nZ, int direction_size, int time, double gamma_dot, double c_s,
                                double* density_field,
                                double* velocity_field,
                                double* previous_particle_distributions,
                                double* particle_distributions,
                                const int* directions,
                                const double* weights
);

static void register_stream_lees_edwards_functions() {
    // TODO johannes: flops
    add_stream_lees_edwards_array_func(&stream_lees_edwards_arrays, "Stream Lees Edwards Baseline - Arrays", 10);
    add_stream_lees_edwards_struct_func(&stream_lees_edwards_baseline, "Stream Lees Edwards Baseline - Structs", 10);
}

#endif //CMDLINE_LBM_LEES_EDWARDS_H
