
#ifndef CMDLINE_LBM_LEES_EDWARDS_H
#define CMDLINE_LBM_LEES_EDWARDS_H
#include "timing.h"

// Flops: 4 + (nX * nZ * q / 3) * 233
// Intops: 1 + nX * nZ * q * (32 * nY + 106 / 3)
// Doubles: // particle_dist: nXYZ*q
            // prev_particle_dist: nXYZ*q;    write
            // velocity_field: nXYZ*3
            // denssity field nXYZ, 
            // weights q
//  Ints:   // dirs: 3xq
void stream_lees_edwards_baseline(struct LBMarrays* S, int time);

void stream_lees_edwards_structural(struct LBMarrays* S, int time);

void stream_lees_edwards_loop_order(struct LBMarrays* S, int time);

void stream_lees_edwards_loop_copy(struct LBMarrays* S, int time);

void stream_lees_edwards_avx(struct LBMarrays* S, int time);

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
    long val0 = S->nX * S->nY * S->nZ;
    struct ops ops = {
            4 + (val / 3 * 233),
            1 + (val * (32 * S->nY + 106 / 3)),
            val0*(S->direction_size*2+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(3*sizeof(int) +sizeof(double)),
            val0*S->direction_size*(int)(sizeof(double))

    };
    return ops;
}


// iops: 1 + i * (9 + 3 * z + 3 * z * y + 13 * z * y * x) + 2 * i/3 * (6 + 6 * z + 24 * z * x)
// flops: 9 + (2 * i/3 * z * x * 71)
// reads: 5 + q * 3 ints, q + 2 + (2 * q + 1)* x * y * z doubles
// writes: q * x * y * z doubles, 0 ints
static struct ops stream_lees_edwards_flops_1(struct LBMarrays* S) {
    long val0 = S->nX * S->nY * S->nZ;
    struct ops ops = {
            9 + (2 * S->direction_size / 3 * S->nZ * S->nX * 71),
            1 + S->direction_size * (9 + 3 * S->nZ + 3 * S->nZ * S->nY + 13 * S->nXYZ) + 2 * S->direction_size * (6 + 6 * S->nZ + 24 * S->nX * S->nZ) / 3,
            (5 + 3 * S->direction_size)*(long)(sizeof(int)) + (S->direction_size + 2 + (2 * S->direction_size + 1) * val0)*(long)(sizeof(double)),
            val0 * S->direction_size*(int)(sizeof(double))
    };
    return ops;
}

static struct ops stream_lees_edwards_flops_2(struct LBMarrays* S) {
    long q = S->direction_size;
    long x = S->nX;
    long y = S->nY;
    long z = S->nZ;
    long grid = S->nX * S->nY * S->nZ;
    struct ops ops = {
            9 + 39 * (q * 2 / 3) * x * z,
            1 + q * (9 + 3 * z + 3 * z * y + 13 * S->nXYZ) + 2 * S->direction_size * (6 + 6 * z + 24 * x * z) / 3,
            (5 + 3 * q)*(long)(sizeof(int)) + (q + 2 + (2 * q + 1) * grid)*(long)(sizeof(double)),
            grid * q *(int)(sizeof(double))
    };
    return ops;
}

static struct ops stream_lees_edwards_flops_3(struct LBMarrays* S) {
    long q = S->direction_size;
    long upDownBoundary = S->direction_size * 2 / 3;
    long upBoundary = S->direction_size / 3;
    long x = S->nX;
    long y = S->nY;
    long z = S->nZ;
    long grid = S->nX * S->nY * S->nZ;


    long periodicIops;
    if(q == 27) {
        periodicIops = 13 + 6 * q + 78 * z + 18 * z * (7 + 12 * (S->nY - 2));
    } else if (q == 15) {
        periodicIops = 13 + 6 * q + 26 * z + 10 * z * (7 + 12 * (S->nY - 2));
    } else {
        periodicIops = 1 + 6 * q + 26 * z +  6 * z * (7 + 12 * (S->nY - 2));
    }

    struct ops ops = {
            9 + 39 * (q * 2 / 3) * x * z,
            8 + upDownBoundary + 3 * upBoundary + 6 * upDownBoundary * z + 22 * upDownBoundary * x * z +  periodicIops,
            (5 + 3 * q)*(long)(sizeof(int)) + (q + 2 + (2 * q + 1) * grid)*(long)(sizeof(double)),
            grid * q *(int)(sizeof(double))
    };
    return ops;
}
static void register_stream_lees_edwards_functions() {
    add_stream_lees_edwards_struct_func(&stream_lees_edwards_structural, &stream_lees_edwards_flops_1, "Stream Lees Edwards - Structual Optimizations");
    add_stream_lees_edwards_struct_func(&stream_lees_edwards_loop_order, &stream_lees_edwards_flops_2, "Stream Lees Edwards - Loop Order");
    add_stream_lees_edwards_struct_func(&stream_lees_edwards_loop_copy, &stream_lees_edwards_flops_3, "Stream Lees Edwards - Memcpy");
    add_stream_lees_edwards_struct_func(&stream_lees_edwards_avx, &stream_lees_edwards_flops_3, "Stream Lees Edwards - AVX");
}

#endif //CMDLINE_LBM_LEES_EDWARDS_H
