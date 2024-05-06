#include "timing.h"

#ifndef CMDLINE_LBM_PERIODIC_H
#define CMDLINE_LBM_PERIODIC_H

void stream_periodic_baseline(struct LBMarrays* S, int time);

void stream_periodic_arrays(int nX, int nY, int nZ, int direction_size,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions
);

static void register_stream_periodic_functions() {
    // TODO johannes: flops
    add_stream_periodic_array_func(&stream_periodic_arrays, "Stream Periodic Baseline - Arrays", 10);
    add_stream_periodic_struct_func(&stream_periodic_baseline, "Stream Periodic Baseline - Structs", 10);
}

#endif //CMDLINE_LBM_PERIODIC_H
