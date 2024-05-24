#include "LBM.h"
#include "string.h"

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
                    int zmd = (S->nZ + z - directionY) % S->nZ;
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
                    double factorPos = weight * S->density_field[shiftedIndex];

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
                    double velocityShiftX1 =velocityShiftX - u_le_x; // 1 Flop
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
            int yDistIndex = distIndex + (S->nY - 1)* S->nX;
            for (int x = 0; x < S->nX; x++) {
                int x_pos = (x - d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                int x_shifted = (x - d_x_I + xDimMinusDirection - 1) % S->nX; // 6 Intops
                for (int z = 0; z < S->nZ; z++) {
                    int zmd = (S->nZ + z - directionY) % S->nZ;
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
                    int zmd = (S->nZ + z - directionY) % S->nZ;
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
    // 9 flops
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
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int yzIndex = ymd * S->nX + zmd * S->nXY;
                for (int x = 0; x < S->nX; x++) {
                    int index = x + z * S->nXY + distIndex;
                    // Bottom Wall.
                    // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                    // by Eric Lorenz and Alfons G. Hoekstra
                    int x_pos = (x + d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                    int x_shifted = (x + d_x_I + xDimMinusDirection + 1) % S->nX; // 6 Intops
                    int posIndex = x_pos + yzIndex; // 5 Intops
                    int shiftedIndex = x_shifted + yzIndex; // 5 Intops

                    double velocityPosY = S->velocity_field[3 * posIndex + 1];
                    double velocityPosZ = S->velocity_field[3 * posIndex + 2];
                    double velocityPosX = S->velocity_field[3 * posIndex]; // 1 Flop
                    double normPosYZ = velocityPosY * velocityPosY + velocityPosZ * velocityPosZ;
                    double dotPosYZ = velocityPosY * directionY + velocityPosZ * directionZ;
                    double factorPos = weight * S->density_field[shiftedIndex];

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
                    double velocityShiftX1 =velocityShiftX - u_le_x; // 1 Flop
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
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int yzIndex = zmd * S->nXY;
                for (int x = 0; x < S->nX; x++) {
                    int index = x + z * S->nXY + yDistIndex;
                    // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                    // by Eric Lorenz and Alfons G. Hoekstra
                    // Top Wall.
                    int x_pos = (x - d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                    int x_shifted = (x - d_x_I + xDimMinusDirection - 1) % S->nX; // 6 Intops
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

        for (int z = 0; z < S->nZ; z++) {
            int zmd = (S->nZ + z - directionY) % S->nZ;
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
        int directionX = S->directions[3 * i ];
        int directionY = S->directions[3 * i + 1];
        int directionZ = S->directions[3 * i + 2];
        int distIndex = i * S->nXYZ; // 3 Intops

        int xDimMinusDirection = 2 * S->nX - directionX;
        double weight = S->weights[i];
        int startY = 0;
        int endY = S->nY;
        if(directionY == 1) {
            int ymd = S->nY - 1;
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int yzIndex = ymd * S->nX + zmd * S->nXY;
                for (int x = 0; x < S->nX; x++) {
                    int index = x + z * S->nXY + distIndex;
                    // Bottom Wall.
                    // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                    // by Eric Lorenz and Alfons G. Hoekstra
                    int x_pos = (x + d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                    int x_shifted = (x + d_x_I + xDimMinusDirection + 1) % S->nX; // 6 Intops
                    int posIndex = x_pos + yzIndex; // 5 Intops
                    int shiftedIndex = x_shifted + yzIndex; // 5 Intops

                    double velocityPosX = S->velocity_field[3 * posIndex];
                    double velocityPosY = S->velocity_field[3 * posIndex + 1];
                    double velocityPosZ = S->velocity_field[3 * posIndex + 2];
                    double normPosYZ = velocityPosY * velocityPosY + velocityPosZ * velocityPosZ;
                    double dotPosYZ = velocityPosY * directionY + velocityPosZ * directionZ;
                    double factorPos = weight * S->density_field[shiftedIndex];

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

                    double velocityShiftX = S->velocity_field[3 * shiftedIndex];
                    double velocityShiftY = S->velocity_field[3 * shiftedIndex + 1];
                    double velocityShiftZ = S->velocity_field[3 * shiftedIndex + 2];
                    double normShiftYZ = velocityShiftY * velocityShiftY + velocityShiftZ * velocityShiftZ;
                    double dotShiftYZ = velocityShiftY * directionY + velocityShiftZ * directionZ;
                    double factorShift = weight * S->density_field[shiftedIndex];

                    // calculate_feq_u 3
                    double velocityShiftX1 =velocityShiftX - u_le_x; // 1 Flop
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
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int yzIndex = zmd * S->nXY;
                for (int x = 0; x < S->nX; x++) {
                    int index = x + z * S->nXY + yDistIndex;
                    // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                    // by Eric Lorenz and Alfons G. Hoekstra
                    // Top Wall.
                    int x_pos = (x - d_x_I + xDimMinusDirection) % S->nX; // 5 Intops
                    int x_shifted = (x - d_x_I + xDimMinusDirection - 1) % S->nX; // 6 Intops
                    int posIndex = x_pos + yzIndex; // 5 Intops
                    int shiftedIndex = x_shifted + yzIndex; // 5 Intops

                    double velocityPosX = S->velocity_field[3 * posIndex];
                    double velocityPosY = S->velocity_field[3 * posIndex + 1];
                    double velocityPosZ = S->velocity_field[3 * posIndex + 2];
                    double normPosYZ = velocityPosY * velocityPosY + velocityPosZ * velocityPosZ;
                    double dotPosYZ = velocityPosY * directionY + velocityPosZ * directionZ;
                    double factorPos = weight * S->density_field[posIndex];

                    // calculate_feq_u 1
                    double velocityPosX1 = S->velocity_fieldX[posIndex] + u_le_x; // 1 Flop
                    double dotPos1 = velocityPosX1 * directionX + dotPosYZ; // 5 Flops
                    double normPos1 = velocityPosX1 * velocityPosX1 + normPosYZ; // 5 Flops
                    double feqPos1 = factorPos * (1.0 + dotPos1 / cSqrt + dotPos1 * dotPos1 / cSqrtSqrtTwo - normPos1 / cSqrtTwo); // 16 Flops

                    // calculate_feq_u 2
                    double dotPos2 = velocityPosX * directionX + dotPosYZ; // 5 Flops
                    double normPos2 = velocityPosX * velocityPosX + normPosYZ; // 5 Flops
                    double feqPos2 = factorPos * (1.0 + dotPos2 / cSqrt + dotPos2 * dotPos2 / cSqrtSqrtTwo - normPos2 / cSqrtTwo); // 16 Flops

                    double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] + feqPos1 - feqPos2;

                    double velocityShiftX = S->velocity_field[3 * shiftedIndex];
                    double velocityShiftY = S->velocity_field[3 * shiftedIndex + 1];
                    double velocityShiftZ = S->velocity_field[3 * shiftedIndex + 2];
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
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + S->nX], (sizeof(double)) * (S->nXY - S->nX));
            }
        } else if(directionX == 0 && directionY == 1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                memcpy(&S->particle_distributions[zIndex + S->nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * (S->nXY - S->nX));
            }
        } else if(directionX == -1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = (S->nY + y - directionY) % S->nY;
                    int otherIndex = ymd * S->nX + otherZIndex;
                    int index = y * S->nX + zIndex;
                    memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[otherIndex + 1], (sizeof(double)) * (S->nX - 1));
                    S->particle_distributions[index + S->nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        } else if(directionX == 1) {
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionY) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = (S->nY + y - directionY) % S->nY;
                    int otherIndex = ymd * S->nX + otherZIndex;
                    int index = y * S->nX + zIndex;
                    memcpy(&S->particle_distributions[index + 1], &S->previous_particle_distributions[otherIndex], (sizeof(double)) * (S->nX - 1));
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex + S->nX - 1];
                }
            }
        }
    }
}