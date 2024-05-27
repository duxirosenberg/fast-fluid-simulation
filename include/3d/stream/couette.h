
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
void stream_couette_code_motion(struct LBMarrays* S);
void stream_couette_loop_structure(struct LBMarrays* S);
void stream_couette_memcpy(struct LBMarrays* S);

void stream_couette_flattened(struct LBMarrays* S);
void stream_couette_flattened_avx(struct LBMarrays* S);



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
            (long)(5 * val / 3),
            (long)(val * (39 * S->nY + 32 / 3)),
            val*S->nY*2*(int)(sizeof(double)) + S->direction_size*(int)((3+1)*sizeof(int) +sizeof(double)),
            val*S->nY*(int)(sizeof(double))
    };
    return ops;
}

static struct ops stream_couette_loop_structure_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    int q = S->direction_size;
    double iops_innermost_loop = (2.0/3 + 14.0/3)+(3.0/3 +14.0/3)+(S->nY-2)*7.0;
    struct ops ops = {
        (long)(4 + q *(S->nZ * S->nX / 3)), //flops
        (long)(q * (11 + S->nZ * ( 2 + S->nX * iops_innermost_loop))), //iops
        (long)(q*(S->nXYZ + 1)*sizeof(double) + 4*q*sizeof(int)), //bytes read
        (long)(q*S->nXYZ*sizeof(double)) //bytes written

    };
    return ops;
}

static struct ops stream_couette_flattened_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    int q = S->direction_size;
    int conditionals = (S->nXYZ + 4*S->nZ*S->nY + 5* S->nZ + 1);
    double iops_innermost_loop = (2.0/3 + 10.0/3)+(3.0/3 +10.0/3)+(S->nY-2)*5.0;
    struct ops ops = {
        (long)(4 + q *(S->nXYZ + 2)), //flops
        (long)(q * (12 + S->nZ * ( 2 + S->nX * iops_innermost_loop)+ conditionals )), //iops
        (long)(q*(S->nXYZ + 1)*sizeof(double) + 4*q*sizeof(int)), //bytes read
        (long)(q*S->nXYZ*sizeof(double)) //bytes written

    };
    return ops;
}





static void register_stream_couette_functions() {
    add_stream_couette_struct_func(&stream_couette_baseline, &stream_couette_baseline_flops, "Stream Couette - Structs Bl");
    add_stream_couette_struct_func(&stream_couette_code_motion, &stream_couette_baseline_flops, "Stream Couette - Structs Code Motion");
    add_stream_couette_struct_func(&stream_couette_loop_structure, &stream_couette_loop_structure_flops, "Stream Couette - Structs Loop Structure");
    add_stream_couette_struct_func(&stream_couette_memcpy, &stream_couette_baseline_flops, "Stream Couette - Structs Memcpy");
    
    //add_stream_couette_struct_func(&stream_couette_flattened, &stream_couette_flattened_flops, "Stream Couette - Structs Loop Flattened");
    //add_stream_couette_struct_func(&stream_couette_flattened_avx, &stream_couette_baseline_flops, "Stream Couette - Structs Loop Structure AVX");

    //add_stream_couette_struct_func(&stream_couette_loop_structure_avx, &stream_couette_baseline_flops, "Stream Couette - Structs Loop Structure AVX");
    //add_stream_couette_array_func(&stream_couette_arrays, &stream_couette_baseline_flops, "Stream Couette - Arrays Bl");
    //add_stream_couette_array_func(&stream_couette_arrays_opt1, &stream_couette_baseline_flops, "Stream Couette - Arrays Opt1: Code Motion and Index Pre-computation");
}

#endif //CMDLINE_LBM_COUETTE_H
