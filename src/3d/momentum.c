#include "LBM.h"

/**
 * These functions calculate the macroscopic moments and the equilibrium distributions as described in
 * the steps 1 & 2 of the Time Step Algorithm described in Section 3.3.2 of The Lattice Boltzmann Method book.
 */


// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
void momentum_baseline(struct LBMarrays* S) {
    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                double new_density = 0;
                double u[] = {0, 0, 0};
                for(int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    double dist = S->particle_distributions[index];
                    new_density += dist;
                    u[0] += dist * (double)S->directions[3 * i];
                    u[1] += dist * (double)S->directions[3 * i + 1];
                    u[2] += dist * (double)S->directions[3 * i + 2];
                }
                int index = (z * S->nX * S->nY) + (y * S->nX) + x;
                S->density_field[index] = new_density;
                S->velocity_field[3 * index] = u[0] / new_density;
                S->velocity_field[3 * index + 1] = u[1] / new_density;
                S->velocity_field[3 * index + 2] = u[2] / new_density;
            }
        }
    }
}

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
void momentum_arrays(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions) {
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            for(int z = 0; z < nZ; z++) {
                //Equation 3.1
                double new_density = 0;
                double u[] = {0, 0, 0};
                for(int i = 0; i < direction_size; i++) {
                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    double dist = particle_distributions[index];
                    new_density += dist;
                    u[0] += dist * directions[3 * i];
                    u[1] += dist * directions[3 * i + 1];
                    u[2] += dist * directions[3 * i + 2];
                }
                int index = (z * nX * nY) + (y * nX) + x;
                density_field[index] = new_density;
                velocity_field[3 * index] = u[0] / new_density;
                velocity_field[3 * index + 1] = u[1] / new_density;
                velocity_field[3 * index + 2] = u[2] / new_density;
            }
        }
    }
}
