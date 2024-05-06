#include "LBM.h"

/**
 * These functions simulate the collision as described in step 3 of the Time Step Algorithm described in
 * Section 3.3.2 of The Lattice Boltzmann Method book.
 */

// Flops: 26
// Intops: 10
static double calculate_feq_arrays(int index,
                                   const double c_s,
                                   double* density_field,
                                   double* velocity_field,
                                   int directionX,
                                   int directionY,
                                   int directionZ,
                                   double weight
) {
    double velocityX = velocity_field[3 * index];
    double velocityY = velocity_field[3 * index + 1];
    double velocityZ = velocity_field[3 * index + 2];

    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
    //Equation 3.4 with c_s^2 = 1/3
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
    return weight * density_field[index] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square / (2 * c_s * c_s));
}

// Flops: 26
// Intops: 10
static double calculate_feq_struct( int feqIndex,
                                    struct LBMarrays* S,
                                    int i
)                      {
    double velocityX = S->velocity_field[3 * feqIndex];
    double velocityY = S->velocity_field[3 * feqIndex + 1];
    double velocityZ = S->velocity_field[3 * feqIndex + 2];
    int directionX = S->directions[3 * i];
    int directionY = S->directions[3 * i + 1];
    int directionZ = S->directions[3 * i + 2];
    double weight = S->weights[i];

    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ; // 5 Flops
    // Equation 3.4 with c_s^2 = 1/3
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ; // 5 Flops
    return weight * S->density_field[feqIndex] * (1.0 + dot_product / (S->c_s * S->c_s) + dot_product * dot_product / (2 * S->c_s * S->c_s * S->c_s * S->c_s) - norm_square / (2 * S->c_s * S->c_s)); // 16 Flops
}

// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
void collision_baseline(struct LBMarrays* S) {
    const double tauinv = 1.0 / S->tau;
    const double omtauinv = 1.0 - tauinv;

    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                for(int i = 0; i < S->direction_size; i++) {
                    int feqIndex = (z * S->nX * S->nY) + (y * S->nX) + x;
                    double feq = calculate_feq_struct(feqIndex, S, i); // 10 Intops, 26 Flops

                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    // Equation 3.9
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tauinv * feq;
                }
            }
        }
    }
}

// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
void collision_arrays(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights
) {//Performs the collision step.
    const double tauinv = 1.0 / tau;
    const double omtauinv = 1.0 - tauinv;

    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            for(int z = 0; z < nZ; z++) {
                for(int i = 0; i < direction_size; i++) {
                    int feqIndex = (z * nX * nY) + (y * nX) + x;
                    double feq = calculate_feq_arrays(feqIndex, c_s, density_field, velocity_field, (double) directions[3 * i], (double) directions[3 * i + 1], (double) directions[3 * i + 2], weights[i]);

                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    //Equation 3.9
                    previous_particle_distributions[index] = omtauinv * particle_distributions[index] + tauinv * feq;

                }
            }
        }
    }
}
