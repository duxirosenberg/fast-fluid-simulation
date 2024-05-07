#include "LBM.h"

// Flops: 27
// Intops: 5
static double calculate_feq_u_struct(int index, double u_le_x, struct LBMarrays* S, int directionX, int directionY, int directionZ, double weight) {
    double velocityX = S->velocity_field[3 * index] + u_le_x; // 1 Flop
    double velocityY = S->velocity_field[3 * index + 1];
    double velocityZ = S->velocity_field[3 * index + 2];
    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ; // 5 Flops
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ; // 5 Flops
    // Equation 3.4 with c_s^2 = 1/3
    return weight * S->density_field[index] * (1.0 + dot_product / (S->c_s * S->c_s) + dot_product * dot_product / (2 * S->c_s * S->c_s * S->c_s * S->c_s) - norm_square / (2 * S->c_s * S->c_s)); // 16 Flops
}

// Flops: 27
// Intops: 5
static double calculate_feq_u_array(int index,
                              double u_le_x,
                              const double c_s,
                              double* density_field,
                              double* velocity_field,
                              int directionX,
                              int directionY,
                              int directionZ,
                              double weight
) {
    double velocityX = velocity_field[3 * index] + u_le_x;
    double velocityY = velocity_field[3 * index + 1];
    double velocityZ = velocity_field[3 * index + 2];
    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
    // TODO johannes: Check if this is actually correct, the directions seem like typos imo
    // double norm_square = velocityX * velocityX + velocityY * directionY + velocityZ * directionZ;
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
    //Equation 3.4 with c_s^2 = 1/3
    return weight * density_field[index] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square / (2 * c_s * c_s));
}




// Flops: 4 + (nX * nZ * q / 3) * 233
// Intops: 1 + nX * nZ * q * (32 * nY + 106 / 3)
void stream_lees_edwards_baseline(struct LBMarrays* S, int time) {
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
                                calculate_feq_u_struct(posIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u_struct(posIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
                        double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u_struct(shiftedIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u_struct(shiftedIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
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
                                calculate_feq_u_struct(posIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u_struct(posIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
                        double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u_struct(shiftedIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u_struct(shiftedIndex, 0, S, directionX, directionY, directionZ, S->weights[i]); // 27 + 27 + 2 Flops, 5 + 5 + 1 Intops
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

// Flops: 4 + (nX * nZ * q / 3) * 233
// Intops: 1 + nX * nZ * q * (32 * nY + 106 / 3)
void stream_lees_edwards_arrays(int nX, int nY, int nZ, int direction_size, int time, double gamma_dot, double c_s,
                                            double* density_field,
                                            double* velocity_field,
                                            double* previous_particle_distributions,
                                            double* particle_distributions,
                                            const int* directions,
                                            const double* weights
) {
    double d_x = gamma_dot * (double)nY * (double)time;
    int d_x_I = d_x;
    double d_x_R = d_x - (double)d_x_I;
    d_x_I = d_x_I % nX;
    double s_1 = d_x_R;
    double s_2 = 1.0 - d_x_R;

    for(int x = 0; x < nX; x++) {
        for (int y = 0; y < nY; y++) {
            for (int z = 0; z < nZ; z++) {
                for (int i = 0; i < direction_size; i++) {
                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    int directionX = directions[3 * i];
                    int directionY = directions[3 * i + 1];
                    int directionZ = directions[3 * i + 2];

                    int ymd = (nY + y - directionY) % nY;
                    int zmd = (nZ + z - directionY) % nZ;
                    if(y==0 && directionY == 1) {
                        //Bottom Wall.
                        //Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                        //by Eric Lorenz and Alfons G. Hoekstra
                        int x_pos = (x + d_x_I + 2 * nX - directionX) % nX;
                        int x_shifted = (x + d_x_I + 2 * nX + 1 - directionX) % nX;
                        int posIndex = x_pos + ymd * nX + zmd * nX * nY;
                        int shiftedIndex = x_shifted + ymd * nX + zmd * nX * nY;
                        int distIndex = i * nX * nY * nZ;

                        double u_le_x = -1 * gamma_dot * (double)nY;
                        double galilean_transformation_pos = previous_particle_distributions[posIndex + distIndex] +
                                calculate_feq_u_array(posIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u_array(posIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
                        double galilean_transformation_shift = previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u_array(shiftedIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u_array(shiftedIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
                        //Equation (18) from same paper.
                        particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos;
                    } else if(y==nY-1 && directionY == -1) {
                        //Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                        //by Eric Lorenz and Alfons G. Hoekstra
                        //Top Wall.
                        int x_pos = (x - d_x_I + 2 * nX - directionX) % nX;
                        int x_shifted = (x - d_x_I + 2 * nX - 1 - directionX) % nX;
                        int posIndex = x_pos + ymd * nX + zmd * nX * nY;
                        int shiftedIndex = x_shifted + ymd * nX + zmd * nX * nY;
                        int distIndex = i * nX * nY * nZ;

                        double u_le_x = gamma_dot * (double)nY;
                        double galilean_transformation_pos = previous_particle_distributions[posIndex + distIndex] +
                                calculate_feq_u_array(posIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u_array(posIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
                        double galilean_transformation_shift = previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u_array(shiftedIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u_array(shiftedIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
                        //Equation (18) from same paper.
                        particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos;
                    } else {
                        //Interior points. (Enforce periodicity along x axis);.
                        int xmd = (nX + x - directionX) % nX;
                        int otherIndex = xmd + ymd * nX + zmd * nX * nY + i * nX * nY * nZ;
                        particle_distributions[index] = previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}