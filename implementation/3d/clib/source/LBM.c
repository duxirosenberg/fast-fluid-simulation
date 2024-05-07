#include "LBM.h"
#include "momentum.h"
#include "collision.h"
#include "couette.h"
#include "periodic.h"
#include "lees_edwards.h"


static void stream_baseline(struct LBMarrays* solver, int time) {
    if (solver->boundary_condition == 1) {
        stream_periodic_baseline(solver, time);
    } else if(solver->boundary_condition == 2) {
        stream_couette_baseline(solver);
    } else if(solver->boundary_condition == 3) {
        stream_lees_edwards_baseline(solver, time);
    }
}

void perform_timestep_baseline(struct LBMarrays* solver, int time) {
    // Flops: nX * nY * nZ * q * 29 + 2
    // Intops: nX * nY * nZ * q * 23
    collision_baseline(solver);

    // Depends on chosen boundary condition :)
    stream_baseline(solver, time);

    // Flops: nX * nY * nZ * (3 + 4 * q)
    // Intops: nX * nY * nZ * (10 + 14 * q)
    momentum_baseline(solver);
}

static void stream_arrays(int nX, int nY, int nZ, int direction_size, int time, double gamma_dot, double c_s, int boundary_condition,
                          double* density_field,
                          double* velocity_field,
                          double* previous_particle_distributions,
                          double* particle_distributions,
                          const int* directions,
                          const double* weights,
                          int* reverse_indexes
) {
    if (boundary_condition == 1) {
        stream_periodic_arrays(nX, nY, nZ, direction_size, previous_particle_distributions, particle_distributions, directions);
    } else if(boundary_condition == 2) {
        stream_couette_arrays(nX, nY, nZ, direction_size, c_s, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);
    } else if(boundary_condition == 3) {
        stream_lees_edwards_arrays(nX, nY, nZ, direction_size, time, gamma_dot, c_s, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights);
    }
}

void perform_timestep_array(int nX, int nY, int nZ, int direction_size, int time, double tau, double gamma_dot, double c_s, int boundary_condition,
                            double* density_field,
                            double* velocity_field,
                            double* previous_particle_distributions,
                            double* particle_distributions,
                            const int* directions,
                            const double* weights,
                            int* reverse_indexes
) {
    // Flops: nX * nY * nZ * q * 29 + 2
    // Intops: nX * nY * nZ * q * 23
    collision_arrays(nX, nY, nZ, direction_size, tau, c_s, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights);

    // Depends on chosen boundary condition :)
    stream_arrays(nX, nY, nZ, direction_size, time, gamma_dot, c_s, boundary_condition, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);

    // Flops: nX * nY * nZ * (3 + 4 * q)
    // Intops: nX * nY * nZ * (10 + 14 * q)
    momentum_arrays(nX, nY, nZ, direction_size, density_field, velocity_field, particle_distributions, directions);
}