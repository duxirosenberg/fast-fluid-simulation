
#ifndef LBM_H_FUNCTIONS_FILE
#define LBM_H_FUNCTIONS_FILE

#include "timing.h"

void perform_timestep_baseline(struct LBMarrays* S, int time);

void perform_timestep_struct(struct LBMarrays* S, int time);

void perform_timestep_array(int nX, int nY, int nZ, int direction_size, int time, double tau, double gamma_dot, double c_s, int boundary_condition,
                            double* density_field,
                            double* velocity_field,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions,
                            const double* weights,
                            int* reverse_indexes);

static void register_functions() {
    // TODO johannes: flops
    add_array_func(&perform_timestep_array, "Baseline - Arrays", 10);
    add_struct_func(&perform_timestep_struct, "Baseline - Structs", 10);
}

#endif