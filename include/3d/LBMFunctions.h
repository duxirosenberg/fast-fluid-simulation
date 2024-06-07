
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


void perform_timestep_1(struct LBMarrays* S, int time);

void perform_timestep_2(struct LBMarrays* S, int time);

void perform_timestep_3(struct LBMarrays* S, int time);

void perform_timestep_4(struct LBMarrays* S, int time);

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

static struct ops perform_timestep_1_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops collisionOps = collision_flops_2(S);
    struct ops momentumOps = momentum_O1_flops(S);
    struct ops streamOps;

    if(S->boundary_condition == 1) {
        streamOps = stream_periodic_O1_flops(S);
    } else if(S->boundary_condition == 2) {
        streamOps = stream_couette_codemotion_flops(S);
    } else if(S->boundary_condition == 3) {
        streamOps = stream_lees_edwards_flops_1(S);
    }

    struct ops ops = {
            collisionOps.flops + momentumOps.flops + streamOps.flops,
            collisionOps.iops + momentumOps.iops + streamOps.iops,
            collisionOps.bytes_read + momentumOps.bytes_read + streamOps.bytes_read,
            collisionOps.bytes_write + momentumOps.bytes_write + streamOps.bytes_write
    };
    return ops;
}


static struct ops perform_timestep_2_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops collisionOps = collision_flops_nb(S);
    struct ops momentumOps = momentum_O2_flops(S);
    struct ops streamOps;

    if(S->boundary_condition == 1) {
        streamOps = stream_periodic_O2_flops(S);
    } else if(S->boundary_condition == 2) {
        streamOps = stream_couette_loop_structure_flops(S);
    } else if(S->boundary_condition == 3) {
        streamOps = stream_lees_edwards_flops_2(S);
    }

    struct ops ops = {
            collisionOps.flops + momentumOps.flops + streamOps.flops,
            collisionOps.iops + momentumOps.iops + streamOps.iops,
            collisionOps.bytes_read + momentumOps.bytes_read + streamOps.bytes_read,
            collisionOps.bytes_write + momentumOps.bytes_write + streamOps.bytes_write
    };
    return ops;
}


static struct ops perform_timestep_3_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops collisionOps = collision_flops_nb(S);
    struct ops momentumOps = momentum_O2_flops(S);
    struct ops streamOps;

    if(S->boundary_condition == 1) {
        streamOps = stream_periodic_memcpy_flops(S);
    } else if(S->boundary_condition == 2) {
        streamOps = stream_couette_memcpy_flops(S);
    } else if(S->boundary_condition == 3) {
        streamOps = stream_lees_edwards_flops_3(S);
    }

    struct ops ops = {
            collisionOps.flops + momentumOps.flops + streamOps.flops,
            collisionOps.iops + momentumOps.iops + streamOps.iops,
            collisionOps.bytes_read + momentumOps.bytes_read + streamOps.bytes_read,
            collisionOps.bytes_write + momentumOps.bytes_write + streamOps.bytes_write
    };
    return ops;
}


static struct ops perform_timestep_4_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops collisionOps = collision_flops_nb(S);
    struct ops momentumOps = momentum_O2_flops(S);
    struct ops streamOps;

    if(S->boundary_condition == 1) {
        streamOps = stream_periodic_memcpy_flops(S);
    } else if(S->boundary_condition == 2) {
        streamOps = stream_couette_memcpy_flops(S);
    } else if(S->boundary_condition == 3) {
        streamOps = stream_lees_edwards_flops_3(S);
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
    add_lbm_struct_func(&perform_timestep_1, &perform_timestep_1_flops, "LBM 1");
    add_lbm_struct_func(&perform_timestep_2, &perform_timestep_2_flops, "LBM 2");
    add_lbm_struct_func(&perform_timestep_3, &perform_timestep_3_flops, "LBM 3");
    add_lbm_struct_func(&perform_timestep_4, &perform_timestep_4_flops, "LBM 4");
}
#endif