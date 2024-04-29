#include "LBM.h"
#include "LBMFunctions.h"

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
static void compute_density_momentum_moment(struct LBMarrays* S) {
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

// Flops: 26
// Intops: 10
static double calculate_feq( int feqIndex,
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

// Flops: 27
// Intops: 5
static double calculate_feq_u(int index, double u_le_x, struct LBMarrays* S, int directionX, int directionY, int directionZ, double weight) {
    double velocityX = S->velocity_field[3 * index] + u_le_x; // 1 Flop
    double velocityY = S->velocity_field[3 * index + 1];
    double velocityZ = S->velocity_field[3 * index + 2];
    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ; // 5 Flops
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ; // 5 Flops
    // Equation 3.4 with c_s^2 = 1/3
    return weight * S->density_field[index] * (1.0 + dot_product / (S->c_s * S->c_s) + dot_product * dot_product / (2 * S->c_s * S->c_s * S->c_s * S->c_s) - norm_square / (2 * S->c_s * S->c_s)); // 16 Flops
}

// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
static void collision(struct LBMarrays* S) {
    const double tauinv = 1.0 / S->tau;
    const double omtauinv = 1.0 - tauinv;

    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                for(int i = 0; i < S->direction_size; i++) {
                    int feqIndex = (z * S->nX * S->nY) + (y * S->nX) + x;
                    double feq = calculate_feq(feqIndex, S, i); // 10 Intops, 26 Flops

                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    // Equation 3.9
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tauinv * feq;
                }
            }
        }
    }
}

// Flops: 0
// Intops: xN * nY * nZ * q * 32
static void periodic_boundary_condition(struct LBMarrays* S) {
    for(int x = 0; x < S->nX; x++) {
        for (int y = 0; y < S->nY; y++) {
            for (int z = 0; z < S->nZ; z++) {
                for (int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ; // 9 Intops
                    int xmd = (S->nX + x - S->directions[3 * i]) % S->nX; // 4 Intops
                    int ymd = (S->nY + y - S->directions[3 * i + 1]) % S->nY; // 5 Intops
                    int zmd = (S->nZ + z - S->directions[3 * i + 2]) % S->nZ; // 5 Intops
                    int otherIndex = xmd + ymd * S->nX + zmd * S->nX * S->nY + i * S->nX * S->nY * S->nZ; // 9 Intops

                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}

// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
static void couette_boundary_condition(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;
    for(int x = 0; x < S->nX; x++) {
        for (int y = 0; y < S->nY; y++) {
            for (int z = 0; z < S->nZ; z++) {
                for (int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ; // 9 Intops
                    int reverseIndex =  x + y * S->nX + z * S->nX * S->nY + S->reverse_indexes[i] * S->nX * S->nY * S->nZ; // 9 Intops
                    int xDirection = S->directions[3 * i];
                    int yDirection = S->directions[3 * i + 1];
                    int zDirection = S->directions[3 * i + 2];
                    if (y == 0 && yDirection == 1) { // 0 Intops, 0 Flops (nX * nZ * direction_size / 3) times
                        // Bottom Wall.
                        // Equation 5.27 from LBM Principles and Practice.
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == S->nY - 1 && yDirection == -1) { // 0 Intops, 5 Flops (nX * nZ * direction_size / 3) times
                        // Top wall
                        // Equation 5.28 from LBM Principles and Practice.
                        // Coefficients of Equation 5.28 calculated from footnote 17.
                        double u_max = 0.1;
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection * 2 * S->weights[i] / (c_s_square) * u_max;
                    } else { // 16 Intops, 0 Flops
                        // Chapter 13 Taylor-Green periodicity from same book.
                        int xmd = (S->nX + x - xDirection) % S->nX;
                        int ymd = y - yDirection;
                        int zmd = (S->nZ + z - zDirection) % S->nZ;
                        int otherIndex = xmd + ymd * S->nX + zmd * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}

// Flops: 4 + (nX * nZ * q / 3) * 233
// Intops: 1 + nX * nZ * q * (32 * nY + 106 / 3)
static void lees_edwards_boundary_condition(struct LBMarrays* S, int time) {
    double d_x = S->gamma_dot * (double)S->nY * (double)time;
    int d_x_I = d_x;
    double d_x_R = d_x - (double)d_x_I;
    d_x_I = d_x_I % S->nX;
    double s_1 = d_x_R;
    double s_2 = 1.0 - d_x_R;

    for(int x = 0; x < S->nX; x++) {
        for (int y = 0; y < S->nY; y++) {
            for (int z = 0; z < S->nZ; z++) {
                for (int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    int directionX = S->directions[3 * i];
                    int directionY = S->directions[3 * i + 1];
                    int directionZ = S->directions[3 * i + 2];

                    int ymd = (S->nY + y - directionY) % S->nY;
                    int zmd = (S->nZ + z - directionY) % S->nZ;
                    if(y==0 && directionY == 1) {
                        // Bottom Wall.
                        // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                        // by Eric Lorenz and Alfons G. Hoekstra
                        int x_pos = (x + d_x_I + 2 * S->nX - directionX) % S->nX; // 5 Intops
                        int x_shifted = (x + d_x_I + 2 * S->nX + 1 - directionX) % S->nX; // 6 Intops
                        int posIndex = x_pos + ymd * S->nX + zmd * S->nX * S->nY; // 5 Intops
                        int shiftedIndex = x_shifted + ymd * S->nX + zmd * S->nX * S->nY; // 5 Intops
                        int distIndex = i * S->nX * S->nY * S->nZ; // 3 Intops

                        double u_le_x = -1 * S->gamma_dot * (double)S->nY;  // 2 Flops
                        double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] +
                                calculate_feq_u(posIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(posIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
                        double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u(shiftedIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(shiftedIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
                        // Equation (18) from the same paper.
                        S->particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos; // 3 Flops
                    } else if(y==S->nY-1 && directionY == -1) {
                        // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                        // by Eric Lorenz and Alfons G. Hoekstra
                        // Top Wall.
                        int x_pos = (x - d_x_I + 2 * S->nX - directionX) % S->nX; // 5 Intops
                        int x_shifted = (x - d_x_I + 2 * S->nX - 1 - directionX) % S->nX; // 6 Intops
                        int posIndex = x_pos + ymd * S->nX + zmd * S->nX * S->nY; // 5 Intops
                        int shiftedIndex = x_shifted + ymd * S->nX + zmd * S->nX * S->nY; // 5 Intops
                        int distIndex = i * S->nX * S->nY * S->nZ; // 3 Intops

                        double u_le_x = S->gamma_dot * (double)S->nY; // 1 Flop
                        double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] +
                                calculate_feq_u(posIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(posIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
                        double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u(shiftedIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(shiftedIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
                        // Equation (18) from the same paper.
                        S->particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos; // 3 Flops
                    } else {
                        // Interior points. (Enforce periodicity along x-axis).
                        int xmd = (S->nX + x - directionX) % S->nX; // 3 Intops
                        int otherIndex = xmd + ymd * S->nX + zmd * S->nX * S->nY + i * S->nX * S->nY * S->nZ; // 9 Intops
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}

static void stream(struct LBMarrays* S, int time) {
    if (S->boundary_condition == 1) {
        periodic_boundary_condition(S);
    } else if (S->boundary_condition == 2) {
        couette_boundary_condition(S);
    } else if (S->boundary_condition == 3) {
        lees_edwards_boundary_condition(S, time);
    }
}

void perform_timestep_baseline(struct LBMarrays* S, int time) {
    // Flops: nX * nY * nZ * q * 29 + 2
    // Intops: nX * nY * nZ * q * 23
    collision(S);
    // Depends on chosen boundary condition :)
    stream(S, time);
    // Flops: nX * nY * nZ * (3 + 4 * q)
    // Intops: nX * nY * nZ * (10 + 14 * q)
    compute_density_momentum_moment(S);
}
