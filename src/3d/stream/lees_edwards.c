#include "LBM.h"
#include "string.h"
#include "stdlib.h"
#include "immintrin.h"

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
                    int zmd = (S->nZ + z - directionZ) % S->nZ;
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
                    int zmd = (nZ + z - directionZ) % nZ;
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

// Summarize Int operations
// Make direction_size the outermost loop, so we have to fetch directions only once
// inline function calculate_feq_u function
// unrolled loop for applying boundary condition on bottom & top wall
void stream_lees_edwards_structural(struct LBMarrays* S, int time) {
    double d_x = S->gamma_dot * (double)S->nY * (double)time;
    int d_x_I = (int) d_x;
    double d_x_R = d_x - (double)d_x_I;
    d_x_I = d_x_I % S->nX;
    double s_1 = d_x_R;
    double s_2 = 1.0 - d_x_R;
    double u_le_x = S->gamma_dot * (double)S->nY;

    double cSqrt = S->c_s * S->c_s;
    double cSqrtTwo = 2 * cSqrt;
    double cSqrtSqrtTwo = 2 * (cSqrt * cSqrt);

    for (int i = 0; i < S->direction_size; i++) {
        int directionX = S->directions[3 * i];
        int directionY = S->directions[3 * i + 1];
        int directionZ = S->directions[3 * i + 2];
        int distIndex = i * S->nXYZ; // 3 Intops

        int xDimMinusDirection = 2 * S->nX - directionX;
        double weight = S->weights[i];
        int startY = 0;
        int endY = S->nY;
        if(directionY == 1) {
            int ymd = S->nY - 1;
            for (int x = 0; x < S->nX; x++) {
                int x_pos = (x + d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                int x_shifted = (x + d_x_I + xDimMinusDirection + 1) % S->nX; // 6 Intops
                for (int z = 0; z < S->nZ; z++) {
                    int zmd = (S->nZ + z - directionZ) % S->nZ;
                    int yzIndex = ymd * S->nX + zmd * S->nXY;
                    int index = x + z * S->nXY + distIndex;
                    // Bottom Wall.
                    // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                    // by Eric Lorenz and Alfons G. Hoekstra
                    int posIndex = x_pos + yzIndex; // 5 Intops
                    int shiftedIndex = x_shifted + yzIndex; // 5 Intops

                    double velocityPosY = S->velocity_field[3 * posIndex + 1];
                    double velocityPosZ = S->velocity_field[3 * posIndex + 2];
                    double velocityPosX = S->velocity_field[3 * posIndex]; // 1 Flop
                    double normPosYZ = velocityPosY * velocityPosY + velocityPosZ * velocityPosZ;
                    double dotPosYZ = velocityPosY * directionY + velocityPosZ * directionZ;
                    double factorPos = weight * S->density_field[posIndex];

                    // calculate_feq_u 1
                    double velocityPosX1 = velocityPosX - u_le_x; // 1 Flop
                    double dotPos1 = velocityPosX1 * directionX + dotPosYZ; // 5 Flops
                    double normPos1 = velocityPosX1 * velocityPosX1 + normPosYZ; // 5 Flops
                    double feqPos1 = factorPos * (1.0 + dotPos1 / cSqrt + dotPos1 * dotPos1 / cSqrtSqrtTwo - normPos1 / cSqrtTwo); // 16 Flops

                    // calculate_feq_u 2
                    double dotPos2 = velocityPosX * directionX + dotPosYZ; // 5 Flops
                    double normPos2 = velocityPosX * velocityPosX + normPosYZ; // 5 Flops
                    double feqPos2 = factorPos * (1.0 + dotPos2 / cSqrt + dotPos2 * dotPos2 / cSqrtSqrtTwo - normPos2 / cSqrtTwo); // 16 Flops

                    double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] + feqPos1 - feqPos2;

                    double velocityShiftY = S->velocity_field[3 * shiftedIndex + 1];
                    double velocityShiftZ = S->velocity_field[3 * shiftedIndex + 2];
                    double velocityShiftX = S->velocity_field[3 * shiftedIndex];
                    double normShiftYZ = velocityShiftY * velocityShiftY + velocityShiftZ * velocityShiftZ;
                    double dotShiftYZ = velocityShiftY * directionY + velocityShiftZ * directionZ;
                    double factorShift = weight * S->density_field[shiftedIndex];

                    // calculate_feq_u 3
                    double velocityShiftX1 = velocityShiftX - u_le_x; // 1 Flop
                    double dotShift1 = velocityShiftX1 * directionX + dotShiftYZ; // 5 Flops
                    double normShift1 = velocityShiftX1 * velocityShiftX1 + normShiftYZ; // 5 Flops
                    double feqShift1 = factorShift * (1.0 + dotShift1 / cSqrt + dotShift1 * dotShift1 / cSqrtSqrtTwo - normShift1 / cSqrtTwo); // 16 Flops

                    // calculate_feq_u 4
                    double dotShift2 = velocityShiftX * directionX + dotShiftYZ; // 5 Flops
                    double normShift2 = velocityShiftX * velocityShiftX + normShiftYZ; // 5 Flops
                    double feqShift2 = factorShift * (1.0 + dotShift2 / cSqrt + dotShift2 * dotShift2 / cSqrtSqrtTwo - normShift2 / cSqrtTwo); // 16 Flops

                    double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] + feqShift1 - feqShift2;
                    // Equation (18) from the same paper.
                    S->particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos; // 3 Flops
                }
            }
            startY++;
        }
        if(directionY == -1) {
            int yDistIndex = distIndex + (S->nY - 1) * S->nX;
            for (int x = 0; x < S->nX; x++) {
                int x_pos = (x - d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                int x_shifted = (x - d_x_I + xDimMinusDirection - 1) % S->nX; // 6 Intops
                for (int z = 0; z < S->nZ; z++) {
                    int zmd = (S->nZ + z - directionZ) % S->nZ;
                    int yzIndex = zmd * S->nXY;
                    int index = x + z * S->nXY + yDistIndex;
                    // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                    // by Eric Lorenz and Alfons G. Hoekstra
                    // Top Wall.
                    int posIndex = x_pos + yzIndex; // 5 Intops
                    int shiftedIndex = x_shifted + yzIndex; // 5 Intops
                    double velocityPosY = S->velocity_field[3 * posIndex + 1];
                    double velocityPosZ = S->velocity_field[3 * posIndex + 2];
                    double velocityPosX = S->velocity_field[3 * posIndex]; // 1 Flop
                    double normPosYZ = velocityPosY * velocityPosY + velocityPosZ * velocityPosZ;
                    double dotPosYZ = velocityPosY * directionY + velocityPosZ * directionZ;
                    double factorPos = weight * S->density_field[posIndex];

                    // calculate_feq_u 1
                    double velocityPosX1 = S->velocity_field[3 * posIndex] + u_le_x; // 1 Flop
                    double dotPos1 = velocityPosX1 * directionX + dotPosYZ; // 5 Flops
                    double normPos1 = velocityPosX1 * velocityPosX1 + normPosYZ; // 5 Flops
                    double feqPos1 = factorPos * (1.0 + dotPos1 / cSqrt + dotPos1 * dotPos1 / cSqrtSqrtTwo - normPos1 / cSqrtTwo); // 16 Flops

                    // calculate_feq_u 2
                    double dotPos2 = velocityPosX * directionX + dotPosYZ; // 5 Flops
                    double normPos2 = velocityPosX * velocityPosX + normPosYZ; // 5 Flops
                    double feqPos2 = factorPos * (1.0 + dotPos2 / cSqrt + dotPos2 * dotPos2 / cSqrtSqrtTwo - normPos2 / cSqrtTwo); // 16 Flops

                    double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] + feqPos1 - feqPos2;

                    double velocityShiftY = S->velocity_field[3 * shiftedIndex + 1];
                    double velocityShiftZ = S->velocity_field[3 * shiftedIndex + 2];
                    double velocityShiftX = S->velocity_field[3 * shiftedIndex];
                    double normShiftYZ = velocityShiftY * velocityShiftY + velocityShiftZ * velocityShiftZ;
                    double dotShiftYZ = velocityShiftY * directionY + velocityShiftZ * directionZ;
                    double factorShift = weight * S->density_field[shiftedIndex];

                    // calculate_feq_u 3
                    double velocityShiftX1 = velocityShiftX + u_le_x; // 1 Flop
                    double dotShift1 = velocityShiftX1 * directionX + dotShiftYZ; // 5 Flops
                    double normShift1 = velocityShiftX1 * velocityShiftX1 + normShiftYZ; // 5 Flops
                    double feqShift1 = factorShift * (1.0 + dotShift1 / cSqrt + dotShift1 * dotShift1 / cSqrtSqrtTwo - normShift1 / cSqrtTwo); // 16 Flops

                    // calculate_feq_u 4
                    double dotShift2 = velocityShiftX * directionX + dotShiftYZ; // 5 Flops
                    double normShift2 = velocityShiftX * velocityShiftX + normShiftYZ; // 5 Flops
                    double feqShift2 = factorShift * (1.0 + dotShift2 / cSqrt + dotShift2 * dotShift2 / cSqrtSqrtTwo - normShift2 / cSqrtTwo); // 16 Flops

                    double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] + feqShift1 - feqShift2;
                    // Equation (18) from the same paper.
                    S->particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos; // 3 Flops
                }
            }
            endY--;
        }

        for(int x = 0; x < S->nX; x++) {
            int xmd = (S->nX + x - directionX) % S->nX; // 3 Intops
            for (int y = startY; y < endY; y++) {
                int ymd = (S->nY + y - directionY) % S->nY;
                for (int z = 0; z < S->nZ; z++) {
                    int zmd = (S->nZ + z - directionZ) % S->nZ;
                    int yzIndex = ymd * S->nX + zmd * S->nXY;
                    int index = x + y * S->nX + z * S->nXY + distIndex;
                    // Interior points. (Enforce periodicity along x-axis).
                    int otherIndex = xmd + yzIndex + distIndex; // 9 Intops
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}

// iops: 1 + i * (9 + 3 * z + 3 * z * y + 13 * z * y * x) + 2 * i/3 * (6 + 6 * z + 24 * z * x)
// flops: 9 + (2 * i/3 * z * x * 71)
// reads: 5 + q * 3 ints, q + 2 + (2 * q + 1)* x * y * z doubles
// writes: q * x * y * z doubles, 0 ints
// Switched x to be the innermost loop, to match array layout
void stream_lees_edwards_loop_order(struct LBMarrays* S, int time) {
    double d_x = S->gamma_dot * (double)S->nY * (double)time;
    int d_x_I = (int) d_x;
    double d_x_R = d_x - (double)d_x_I;
    d_x_I = d_x_I % S->nX;
    double s_1 = d_x_R;
    double s_2 = 1.0 - d_x_R;
    double u_le_x = S->gamma_dot * (double)S->nY;

    double cSqrt = S->c_s * S->c_s;
    double cSqrtTwo = 2 * cSqrt;
    double cSqrtSqrtTwo = 2 * (cSqrt * cSqrt);
    double* gallileanTransform = malloc(S->nX * sizeof (double));

    for (int i = 0; i < S->direction_size; i++) {
        int directionX = S->directions[3 * i];
        int directionY = S->directions[3 * i + 1];
        int directionZ = S->directions[3 * i + 2];
        int distIndex = i * S->nXYZ; // 3 Intops

        int xDimMinusDirection = 2 * S->nX - directionX;
        double weight = S->weights[i];
        int startY = 0;
        int endY = S->nY;
        if(directionY != 0) {
            int ymd;
            int yDistIndex;
            if(directionY == 1) {
                yDistIndex = distIndex;
                startY++;
                ymd = S->nY - 1;
            } else {
                yDistIndex = distIndex + (S->nY - 1) * S->nX;
                endY--;
                ymd = 0;
            }
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int yzIndex = ymd * S->nX + zmd * S->nXY;

                for(int x = 0; x < S->nX; x++) {
                    int index = x + yzIndex;
                    double velocityX = S->velocity_field[3 * index]; // 1 Flop
                    double velocityY = S->velocity_field[3 * index + 1];
                    double velocityZ = S->velocity_field[3 * index + 2];
                    double normYZ = velocityY * velocityY + velocityZ * velocityZ;
                    double dotYZ = velocityY * directionY + velocityZ * directionZ;
                    double factor = weight * S->density_field[index];
                    double velocityX1 = velocityX - directionY * u_le_x; // 1 Flop
                    double dot1 = velocityX1 * directionX + dotYZ; // 5 Flops
                    double norm1 = velocityX1 * velocityX1 + normYZ; // 5 Flops
                    double feq1 = factor * (1.0 + dot1 / cSqrt + dot1 * dot1 / cSqrtSqrtTwo - norm1 / cSqrtTwo); // 16 Flops
                    double dot2 = velocityX * directionX + dotYZ; // 5 Flops
                    double norm2 = velocityX * velocityX + normYZ; // 5 Flops
                    double feq2 = factor * (1.0 + dot2 / cSqrt + dot2 * dot2 / cSqrtSqrtTwo - norm2 / cSqrtTwo); // 16 Flops

                    gallileanTransform[x] = S->previous_particle_distributions[index + distIndex] + feq1 - feq2;
                }

                for (int x = 0; x < S->nX; x++) {
                    int index = x + z * S->nXY + yDistIndex;
                    int x_pos = (x + directionY * d_x_I + xDimMinusDirection) % S->nX;
                    int x_shifted = (x + directionY * d_x_I + xDimMinusDirection + directionY) % S->nX;
                    S->particle_distributions[index] = s_1 * gallileanTransform[x_shifted] + s_2 * gallileanTransform[x_pos]; // 3 Flops
                }
            }
        }

        for (int z = 0; z < S->nZ; z++) {
            int zmd = (S->nZ + z - directionZ) % S->nZ;
            for (int y = startY; y < endY; y++) {
                int ymd = (S->nY + y - directionY) % S->nY;
                for(int x = 0; x < S->nX; x++) {
                    int yzIndex = ymd * S->nX + zmd * S->nXY;
                    int index = x + y * S->nX + z * S->nXY + distIndex;
                    // Interior points. (Enforce periodicity along x-axis).
                    int xmd = (S->nX + x - directionX) % S->nX; // 3 Intops
                    int otherIndex = xmd + yzIndex + distIndex; // 9 Intops
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}

void stream_lees_edwards_loop_copy(struct LBMarrays* S, int time) {
    // 9 flops
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    int nXY = S->nXY;
    int nXYZ = S->nXYZ;
    double d_x = S->gamma_dot * (double)nY * (double)time;
    int d_x_I = (int) d_x;
    double d_x_R = d_x - (double)d_x_I;
    d_x_I = d_x_I % nX;
    double s_1 = d_x_R;
    double s_2 = 1.0 - d_x_R;
    double u_le_x = S->gamma_dot * (double)nY;

    double cSqrt = S->c_s * S->c_s;
    double cSqrtTwo = 2 * cSqrt;
    double cSqrtSqrtTwo = 2 * (cSqrt * cSqrt);

    double* gallileanTransform = malloc(nX * sizeof (double));


    for (int i = 0; i < S->direction_size; i++) {
        int xDirection = S->directions[3 * i];
        int yDirection = S->directions[3 * i + 1];
        int zDirection = S->directions[3 * i + 2];
        int dirIndex = i * nXYZ; // 3 Intops

        int xDimMinusDirection = 2 * nX - xDirection;
        double weight = S->weights[i];
        int startY = 0;
        int endY = nY;
        if (yDirection != 0) {
            int ymd;
            int yDistIndex;
            if (yDirection == 1) {
                yDistIndex = dirIndex;
                startY++;
                ymd = nY - 1;
            } else {
                yDistIndex = dirIndex + (nY - 1) * nX;
                endY--;
                ymd = 0;
            }
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int yzIndex = ymd * nX + zmd * nXY;

                for (int x = 0; x < nX; x++) {
                    int index = x + yzIndex;
                    double velocityX = S->velocity_fieldX[index]; // 1 Flop
                    double velocityY = S->velocity_fieldY[index];
                    double velocityZ = S->velocity_fieldZ[index];
                    double normYZ = velocityY * velocityY + velocityZ * velocityZ;
                    double dotYZ = velocityY * yDirection + velocityZ * zDirection;
                    double factor = weight * S->density_field[index];
                    double velocityX1 = velocityX - yDirection * u_le_x; // 1 Flop
                    double dot1 = velocityX1 * xDirection + dotYZ; // 5 Flops
                    double norm1 = velocityX1 * velocityX1 + normYZ; // 5 Flops
                    double feq1 =
                            factor * (1.0 + dot1 / cSqrt + dot1 * dot1 / cSqrtSqrtTwo - norm1 / cSqrtTwo); // 16 Flops
                    double dot2 = velocityX * xDirection + dotYZ; // 5 Flops
                    double norm2 = velocityX * velocityX + normYZ; // 5 Flops
                    double feq2 =
                            factor * (1.0 + dot2 / cSqrt + dot2 * dot2 / cSqrtSqrtTwo - norm2 / cSqrtTwo); // 16 Flops

                    gallileanTransform[x] = S->previous_particle_distributions[index + dirIndex] + feq1 - feq2;
                }

                for (int x = 0; x < nX; x++) {
                    int index = x + z * nXY + yDistIndex;
                    int x_pos = (x + yDirection * d_x_I + xDimMinusDirection) % nX;
                    int x_shifted = (x + yDirection * d_x_I + xDimMinusDirection + yDirection) % nX;
                    S->particle_distributions[index] =
                            s_1 * gallileanTransform[x_shifted] + s_2 * gallileanTransform[x_pos]; // 3 Flops
                }
            }
        }

        if (xDirection == 0 && yDirection == 0 && zDirection == 0) {
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex],
                   (sizeof(double)) * nXYZ);
        } else if (xDirection == 0 && yDirection == 0 && zDirection == -1) {
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex + nXY],
                   (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dirIndex + nXYZ - nXY],
                   &S->previous_particle_distributions[dirIndex], (sizeof(double)) * nXY);
        } else if (xDirection == 0 && yDirection == 0 && zDirection == 1) {
            memcpy(&S->particle_distributions[dirIndex + nXY], &S->previous_particle_distributions[dirIndex],
                   (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dirIndex],
                   &S->previous_particle_distributions[dirIndex + nXYZ - nXY], (sizeof(double)) * nXY);
        } else if (xDirection == 0 && yDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + nX],
                       (sizeof(double)) * (nXY - nX));
            }
        } else if (xDirection == 0 && yDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                memcpy(&S->particle_distributions[zIndex + nX], &S->previous_particle_distributions[otherZIndex],
                       (sizeof(double)) * (nXY - nX));
            }
        } else if (xDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = (nY + y - yDirection) % nY;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[otherIndex + 1],
                           (sizeof(double)) * (nX - 1));
                    S->particle_distributions[index + nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        } else if (xDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = (nY + y - yDirection) % nY;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    memcpy(&S->particle_distributions[index + 1], &S->previous_particle_distributions[otherIndex],
                           (sizeof(double)) * (nX - 1));
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex + nX - 1];
                }
            }
        }
    }
}

void stream_lees_edwards_avx(struct LBMarrays* S, int time) {
    double d_x = S->gamma_dot * (double)S->nY * (double)time;
    int d_x_I = (int) d_x;
    double d_x_R = d_x - (double)d_x_I;
    d_x_I = d_x_I % S->nX;
    double s_1 = d_x_R;
    double s_2 = 1.0 - d_x_R;
    double u_le_x = S->gamma_dot * (double)S->nY;

    double cSqrt = S->c_s * S->c_s;
    double cSqrtTwo = 2 * cSqrt;
    double cSqrtSqrtTwo = 2 * (cSqrt * cSqrt);

    __m256d vS1 = _mm256_set1_pd(s_1);
    __m256d vS2 = _mm256_set1_pd(s_2);
    __m256d vOne = _mm256_set1_pd(1.0);
    __m256d vUleX = _mm256_set1_pd(u_le_x);
    __m256d vCSqrt = _mm256_div_pd(vOne, _mm256_set1_pd(cSqrt));
    __m256d vCSqrtTwo = _mm256_div_pd(vOne, _mm256_set1_pd(cSqrtTwo));
    __m256d vCSqrtSqrtTwo = _mm256_div_pd(vOne, _mm256_set1_pd(cSqrtSqrtTwo));

    double* gallileanTransform = malloc(2 * (S->nX + 4) * sizeof (double));

    for (int i = 0; i < S->direction_size; i++) {
        int directionX = S->directions[3 * i ];
        int directionY = S->directions[3 * i + 1];
        int directionZ = S->directions[3 * i + 2];
        int distIndex = i * S->nXYZ; // 3 Intops

        int startY = 0;
        int endY = S->nY;
        if(directionY != 0) {
            double weight = S->weights[i];
            int xDimMinusDirection = 2 * S->nX - directionX;
            __m256d vWeight = _mm256_set1_pd(weight);
            __m256d vDirectionX = _mm256_set1_pd(directionX);
            __m256d vDirectionY = _mm256_set1_pd(directionY);
            __m256d vDirectionZ = _mm256_set1_pd(directionZ);
            int ymd;
            int yDistIndex;
            if(directionY == 1) {
                yDistIndex = distIndex;
                startY++;
                ymd = S->nY - 1;
            } else {
                yDistIndex = distIndex + (S->nY - 1) * S->nX;
                endY--;
                ymd = 0;
            }
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int yzIndex = ymd * S->nX + zmd * S->nXY;

                for(int x = 0; x < S->nX; x += 4) {
                    int index = x + yzIndex;
                    __m256d vY = _mm256_loadu_pd(&S->velocity_fieldY[index]);
                    __m256d vZ = _mm256_loadu_pd(&S->velocity_fieldZ[index]);
                    __m256d normXY = _mm256_fmadd_pd(vZ, vZ, _mm256_mul_pd(vY, vY));
                    __m256d dotXY = _mm256_fmadd_pd(vZ, vDirectionZ, _mm256_mul_pd(vY, vDirectionY));

                    __m256d vX = _mm256_loadu_pd(&S->velocity_fieldX[index]);
                    __m256d vX1 = _mm256_sub_pd(vX, _mm256_mul_pd(vDirectionY, vUleX));// _mm256_fnmsub_pd(vDirectionY, vUleX, vX);
                    __m256d norm1 = _mm256_fmadd_pd(vX1, vX1, normXY);
                    __m256d dot1 = _mm256_fmadd_pd(vX1, vDirectionX, dotXY);
                    __m256d vDensity = _mm256_loadu_pd(&S->density_field[index]);
                    __m256d vFactor = _mm256_mul_pd(vWeight, vDensity);
                    __m256d feq1 = _mm256_mul_pd(
                            vFactor,
                            _mm256_sub_pd(
                                    _mm256_add_pd(vOne, _mm256_fmadd_pd(dot1, _mm256_mul_pd(dot1, vCSqrtSqrtTwo), _mm256_mul_pd(dot1, vCSqrt))),
                                    _mm256_mul_pd(norm1, vCSqrtTwo))
                    );
                    __m256d norm2 = _mm256_fmadd_pd(vX, vX, normXY);
                    __m256d dot2 = _mm256_fmadd_pd(vX, vDirectionX, dotXY);
                    __m256d feq2 = _mm256_mul_pd(
                            vFactor,
                            _mm256_sub_pd(
                                    _mm256_add_pd(vOne, _mm256_fmadd_pd(dot2, _mm256_mul_pd(dot2, vCSqrtSqrtTwo), _mm256_mul_pd(dot2, vCSqrt))),
                                    _mm256_mul_pd(norm2, vCSqrtTwo))
                            );
                    __m256d galileanTransV = _mm256_sub_pd(_mm256_add_pd(_mm256_loadu_pd(&S->previous_particle_distributions[index + distIndex]), feq1), feq2);
                    _mm256_storeu_pd(&gallileanTransform[x], galileanTransV);
                }
                _mm256_storeu_pd(&gallileanTransform[S->nX], _mm256_loadu_pd(&gallileanTransform[0]));
                int index = z * S->nXY + yDistIndex;
                int target = S->nX + z * S->nXY + yDistIndex;
                int x_pos = (directionY * d_x_I + xDimMinusDirection);
                int x_shifted = (directionY * d_x_I + xDimMinusDirection + directionY);
                for (; index < target; index+=4, x_pos+=4, x_shifted+=4) {
                    __m256d transform1 = _mm256_loadu_pd(&gallileanTransform[x_shifted % S->nX]);
                    __m256d transform2 = _mm256_loadu_pd(&gallileanTransform[x_pos % S->nX]);
                    __m256d val = _mm256_fmadd_pd(vS1, transform1, _mm256_mul_pd(vS2, transform2));
                    _mm256_storeu_pd(&S->particle_distributions[index], val);
                }
            }
        }

        if(directionX == 0 && directionY == 0 && directionZ == 0) {
            memcpy(&S->particle_distributions[distIndex], &S->previous_particle_distributions[distIndex], (sizeof(double)) * S->nXYZ);
        } else if(directionX == 0 && directionY == 0 && directionZ == -1) {
            memcpy(&S->particle_distributions[distIndex], &S->previous_particle_distributions[distIndex + S->nXY], (sizeof(double)) * (S->nXYZ - S->nXY));
            memcpy(&S->particle_distributions[distIndex + S->nXYZ - S->nXY], &S->previous_particle_distributions[distIndex], (sizeof(double)) * S->nXY);
        } else if(directionX == 0 && directionY == 0 && directionZ == 1) {
            memcpy(&S->particle_distributions[distIndex + S->nXY], &S->previous_particle_distributions[distIndex], (sizeof(double)) * (S->nXYZ - S->nXY));
            memcpy(&S->particle_distributions[distIndex], &S->previous_particle_distributions[distIndex + S->nXYZ - S->nXY], (sizeof(double)) * S->nXY);
        } else if(directionX == 0 && directionY == -1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + S->nX], (sizeof(double)) * (S->nXY - S->nX));
            }
        } else if(directionX == 0 && directionY == 1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                memcpy(&S->particle_distributions[zIndex + S->nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * (S->nXY - S->nX));
            }
        } else if(directionX == -1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = y - directionY;
                    int otherIndex = ymd * S->nX + otherZIndex;
                    int index = y * S->nX + zIndex;
                    memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[otherIndex + 1], (sizeof(double)) * (S->nX - 1));
                    S->particle_distributions[index + S->nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        } else if(directionX == 1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = y - directionY;
                    int otherIndex = ymd * S->nX + otherZIndex;
                    int index = y * S->nX + zIndex;
                    memcpy(&S->particle_distributions[index + 1], &S->previous_particle_distributions[otherIndex], (sizeof(double)) * (S->nX - 1));
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex + S->nX - 1];
                }
            }
        }
    }
    free(gallileanTransform);
}