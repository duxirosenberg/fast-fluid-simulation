#ifndef CMDLINE_LBM_TIMING_UTILS_H
#define CMDLINE_LBM_TIMING_UTILS_H

#include "LBM.h"

struct ops {
    long flops;
    long iops;
};

typedef void(*comp_func_struct_time)(struct LBMarrays*, int);

typedef void(*comp_func_struct)(struct LBMarrays*);

typedef void(*comp_momentum_arrays)(int, int, int, int, double*, double*, double*, const int*);

typedef void(*comp_collision_arrays)(int, int, int, int, double, double, double*, double*, double*, double*, const int*, const double*);

typedef void(*comp_stream_periodic_arrays)(int, int, int, int, double*, double*, const int*);

typedef void(*comp_stream_couette_arrays)(int, int, int, int, double, double*, double*, const int*, const double*, const int*);

typedef void(*comp_stream_lees_edwards_arrays)(int, int, int, int, int, double, double, double*, double*, double*, double*, const int*, const double*);

typedef void(*comp_func_arrays)(int, int, int, int, int, double, double, double, int, double*, double*, double*, double*, const int*, const double*, int*);


typedef struct ops(*calc_flops)(struct LBMarrays*);

// Function prototypes for registering and performing simulations
void add_momentum_array_func(comp_momentum_arrays f, calc_flops ff, const char* name);
void add_collision_array_func(comp_collision_arrays f, calc_flops ff, const char* name);
void add_stream_periodic_array_func(comp_stream_periodic_arrays f, calc_flops ff, const char* name);
void add_stream_couette_array_func(comp_stream_couette_arrays f, calc_flops ff, const char* name);
void add_stream_lees_edwards_array_func(comp_stream_lees_edwards_arrays f, calc_flops ff, const char* name);
void add_lbm_array_func(comp_func_arrays f, calc_flops ff,  const char* name);

void add_momentum_struct_func(comp_func_struct f, calc_flops ff, const char* name);
void add_collision_struct_func(comp_func_struct f, calc_flops ff, const char* name);
void add_stream_periodic_struct_func(comp_func_struct_time f, calc_flops ff, const char* name);
void add_stream_couette_struct_func(comp_func_struct f, calc_flops ff, const char* name);
void add_stream_lees_edwards_struct_func(comp_func_struct_time f, calc_flops ff, const char* name);
void add_lbm_struct_func(comp_func_struct_time f, calc_flops ff, const char* name);


#endif //CMDLINE_LBM_TIMING_UTILS_H
