#ifndef CMDLINE_LBM_TIMING_UTILS_H
#define CMDLINE_LBM_TIMING_UTILS_H

#include "LBM.h"

typedef void(*comp_func_struct)(struct LBMarrays*, int);

typedef void(*comp_func_arrays)(int, int, int, int, int, double, double, double, int, double*, double*, double*, double*, const int*, const double*, int*);

// Function prototypes for registering and performing simulations
void add_array_func(comp_func_arrays f, const char* name, int flops);
void add_struct_func(comp_func_struct f, const char* name, int flops);


#endif //CMDLINE_LBM_TIMING_UTILS_H
