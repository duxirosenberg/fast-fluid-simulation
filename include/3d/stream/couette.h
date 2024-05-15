
#ifndef CMDLINE_LBM_COUETTE_H
#define CMDLINE_LBM_COUETTE_H
#include "timing.h"

// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
// Doubles: // particle_dist: nXYZ*q        write
            // prev_particle_dist: nXYZ*q;
            // weights q
//  Ints:   // dirs: 3xq
            // rev_indexes: q
void stream_couette_baseline(struct LBMarrays* S);

void stream_couette_opt1(struct LBMarrays* S);
void stream_couette_opt2(struct LBMarrays* S);
//void stream_couette_opt3(struct LBMarrays* S);
//void stream_couette_opt4(struct LBMarrays* S);

void stream_couette_arrays(int nX, int nY, int nZ, int direction_size, double c_s,
                           double* previous_particle_distributions,
                           double* particle_distributions,
                           const int* directions,
                           const double* weights,
                           const int* reverse_indexes);

void stream_couette_arrays_opt1(int nX, int nY, int nZ, int direction_size, double c_s,
                           double* previous_particle_distributions,
                           double* particle_distributions,
                           const int* directions,
                           const double* weights,
                           const int* reverse_indexes);

static struct ops stream_couette_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    struct ops ops = {
            5 * val / 3,
            val * (39 * S->nY + 32 / 3),
            val*S->nY*2*(int)(sizeof(double)) + S->direction_size*(int)((3+1)*sizeof(int) +sizeof(double)),
            val*S->nY*(int)(sizeof(double))

    };
    return ops;
}

static void register_stream_couette_functions() {
    //add_stream_couette_struct_func(&stream_couette_baseline, &stream_couette_baseline_flops, "Stream Couette - Structs Bl");
    add_stream_couette_struct_func(&stream_couette_opt1, &stream_couette_baseline_flops, "Stream Couette - Structs Opt1: Code Motion and Index Pre-computation");
    add_stream_couette_struct_func(&stream_couette_opt2, &stream_couette_baseline_flops, "Stream Couette - Structs Opt2: Loop Order / Structure ");
    //add_stream_couette_struct_func(&stream_couette_opt3, &stream_couette_baseline_flops, "Stream Couette - Structs Opt3: flat loop unoptimized");
    //add_stream_couette_struct_func(&stream_couette_opt4, &stream_couette_baseline_flops, "Stream Couette - Structs Opt3: flat loop optimized");

    //add_stream_couette_array_func(&stream_couette_arrays, &stream_couette_baseline_flops, "Stream Couette - Arrays Bl");
    //add_stream_couette_array_func(&stream_couette_arrays_opt1, &stream_couette_baseline_flops, "Stream Couette - Arrays Opt1: Code Motion and Index Pre-computation");
}

#endif //CMDLINE_LBM_COUETTE_H
