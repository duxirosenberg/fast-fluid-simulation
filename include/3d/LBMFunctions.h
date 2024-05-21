
#ifndef LBM_H_FUNCTIONS_FILE
#define LBM_H_FUNCTIONS_FILE

#include "couette.h"
#include "lees_edwards.h"
#include "periodic.h"
#include "collision.h"
#include "momentum.h"
#include "timing.h"

void perform_timestep_baseline(struct LBMarrays* S, int time);

void perform_timestep_array(int nX, int nY, int nZ, int direction_size, int time, double tau, double gamma_dot, double c_s, int boundary_condition,
                            double* density_field,
                            double* velocity_field,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions,
                            const double* weights,
                            int* reverse_indexes);

static struct ops perform_timestep_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops collisionOps = collision_baseline_flops(S);
    struct ops momentumOps = momentum_baseline_flops(S);
    struct ops streamOps;

    if(S->boundary_condition == 1) {
        streamOps = stream_periodic_baseline_flops(S);
    } else if(S->boundary_condition == 2) {
        streamOps = stream_couette_baseline_flops(S);
    } else if(S->boundary_condition == 3) {
        streamOps = stream_lees_edwards_baseline_flops(S);
    }

    struct ops ops = {
            collisionOps.flops + momentumOps.flops + streamOps.flops,
            collisionOps.iops + momentumOps.iops + streamOps.iops,
            collisionOps.bytes_read + momentumOps.bytes_read + streamOps.bytes_read,
            collisionOps.bytes_write + momentumOps.bytes_write + streamOps.bytes_write
    };
    return ops;
}

static void register_lbm_function() {
    // add_lbm_array_func(&perform_timestep_array, &perform_timestep_baseline_flops, "LBM - Arrays Bl");
    add_lbm_struct_func(&perform_timestep_baseline, &perform_timestep_baseline_flops, "LBM - Structs Bl");
}
#endif