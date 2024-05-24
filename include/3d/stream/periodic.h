
#ifndef CMDLINE_LBM_PERIODIC_H
#define CMDLINE_LBM_PERIODIC_H
#include "timing.h"

// Flops: 0
// Intops: xN * nY * nZ * q * 32
// Doubles: // particle_dist: nXYZ*q;  write
            // pev_particle_dist: nXYZ*q;
// Ints     // dirs: 3xq
void stream_periodic_baseline(struct LBMarrays* S, int time);
void stream_periodic_O1(struct LBMarrays* S, int time);
void stream_periodic_O2(struct LBMarrays* S, int time);
void stream_periodic_O3(struct LBMarrays* S, int time);
void stream_periodic_memcpy(struct LBMarrays* S, int time);
void stream_periodic_arrays(int nX, int nY, int nZ, int direction_size,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions
);


//void stream_periodic_ZYXI(struct LBMarrays* S, int time);
//void stream_periodic_O21(struct LBMarrays* S, int time);
//void stream_periodic_loopunrolling(struct LBMarrays* S, int time);

static struct ops stream_periodic_baseline_flops(struct LBMarrays* S) {

    long val = S->nX * S->nY * S->nZ* S->direction_size;
    struct ops ops = {
            0,
            val * 32,
            val*2*(int)(sizeof(double)) + 3*S->direction_size*(int)(sizeof(int)),
            val*(int)(sizeof(double))
    };
    return ops;
}
// Flops: 0
// Intops: nZ(1 + nY(1 + (nX * ( 2 + 19*q))))
static struct ops stream_periodic_O1_flops(struct LBMarrays* S) {

    long val = S->nX * S->nY * S->nZ* S->direction_size;
    struct ops ops = {
            0,
            S->nZ + S->nZ*S->nY*(1 + S->nX*(2 + 19 * S->direction_size)),
            val*2*(int)(sizeof(double)) + 3*S->direction_size*(int)(sizeof(int)),
            val*(int)(sizeof(double))
    };
    return ops;
}

// Flops: 0
// Intops: 2 + q * (4 + (nZ * (4 + nY * (4 + 6 * nX))))

static struct ops stream_periodic_O2_flops(struct LBMarrays* S) {

    long val = S->nX * S->nY * S->nZ* S->direction_size;
    struct ops ops = {
            0,
            2 + S->direction_size*(4 + S->nZ *(4 + S->nY*(4 + 6 *S->nX))),
            val*2*(int)(sizeof(double)) + 3*S->direction_size*(int)(sizeof(int)),
            val*(int)(sizeof(double))
    };
    return ops;
}

static void register_stream_periodic_functions() {
    add_stream_periodic_struct_func(&stream_periodic_O1, &stream_periodic_O1_flops, "Stream Periodic - O1 - precalc");
    add_stream_periodic_struct_func(&stream_periodic_O2, &stream_periodic_O2_flops, "Stream Periodic - O2 - loop reordering");
    add_stream_periodic_struct_func(&stream_periodic_memcpy, &stream_periodic_O2_flops, "Stream Periodic - Memcopy");
    //add_stream_periodic_struct_func(&stream_periodic_O3, &stream_periodic_baseline_flops, "Stream Periodic - O3 - AVX");
    //add_stream_periodic_struct_func(&stream_periodic_ZYXI, &stream_periodic_baseline_flops, "Stream Periodic - ZYXI loop order");
    //add_stream_periodic_struct_func(&stream_periodic_O21, &stream_periodic_baseline_flops, "Stream Periodic - O21");
    //add_stream_periodic_struct_func(&stream_periodic_loopunrolling, &stream_periodic_baseline_flops, "Stream Periodic - O21");
}

#endif //CMDLINE_LBM_PERIODIC_H
