#include "LBM.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


#ifdef LBM_STRUCT

void initialize(int nX, int nY, int nZ, int direction_size,
                double* density_field,
                double* velocity_field,
                double* previous_particle_distributions,
                double* particle_distributions,
                const int* directions,
                const double* weights,
                int* reverse_indexes
                ) {


    for(int i = 0; i < direction_size; i++) {
        for(int j = 0; j < direction_size; j++) {
            if(directions[3 * i] == -directions[3 * j] && directions[3 * i + 1] == -directions[3 * j + 1] && directions[3 * i + 2] == -directions[3 * j + 2]) {
                reverse_indexes[i] = j;
            }
        }
    }

    int box_flatten_length = nX * nY * nZ;

    for(int i = 0; i < box_flatten_length; i++) {
        density_field[i] = 1;
    }

    for(int i = 0; i < 3 * box_flatten_length; i++) {
        velocity_field[i] = 0;
    }

    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            for(int z = 0; z < nZ; z++) {
                for(int i = 0; i < direction_size; i++) {
                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    previous_particle_distributions[index] = weights[i];
                    particle_distributions[index] = weights[i];
                }
            }
        }
    }
}
// void initializeS(struct LBM* S) {
//     for(int i = 0; i < S->direction_size; i++) {
//         for(int j = 0; j < S->direction_size; j++) {
//             if(S->directions[3 * i] == -S->directions[3 * j] && S->directions[3 * i + 1] == -S->directions[3 * j + 1] && S->directions[3 * i + 2] == -S->directions[3 * j + 2]) {
//                 S->reverse_indexes[i] = j;
//             }
//         }
//     }

//     int box_flatten_length = S->nX * S->nY * S->nZ;

//     for(int i = 0; i < box_flatten_length; i++) {
//         S->density_field[i] = 1;
//     }

//     for(int i = 0; i < 3 * box_flatten_length; i++) {
//         S->velocity_field[i] = 0;
//     }

//     for(int x = 0; x < S->nX; x++) {
//         for(int y = 0; y < S->nY; y++) {
//             for(int z = 0; z < S->nZ; z++) {
//                 for(int i = 0; i < S->direction_size; i++) {
//                     int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
//                     S->previous_particle_distributions[index] = S->weights[i];
//                     S->particle_distributions[index] = S->weights[i];
//                 }
//             }
//         }
//     }
// }

void compute_density_momentum_moment(struct LBM* S) {
    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                // printf("a\n");
                double new_density = 0;
                double u[] = {0, 0, 0};
                for(int i = 0; i < S->direction_size; i++) {
                    // printf("%i a\n", i);
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

double calculate_feq( int feqIndex,
                      struct LBM* S,
                      int i
                     ) {
    double velocityX = S->velocity_field[3 * feqIndex];
    double velocityY = S->velocity_field[3 * feqIndex + 1];
    double velocityZ = S->velocity_field[3 * feqIndex + 2];
    int directionX = (double)(S->directions[3 * i]);
    int directionY = (double)(S->directions[3 * i + 1]);
    int directionZ = (double)(S->directions[3 * i + 2]); 
    double weight = S->weights[i];

    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
    // Equation 3.4 with c_s^2 = 1/3
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
    return weight * S->density_field[feqIndex] * (1.0 + dot_product / (S->c_s * S->c_s) + dot_product * dot_product / (2 * S->c_s * S->c_s * S->c_s * S->c_s) - norm_square / (2 * S->c_s * S->c_s));
}


double calculate_feq_u(int index, double u_le_x, struct LBM* S, int directionX, int directionY, int directionZ, double weight) {
    double velocityX = S->velocity_field[3 * index] + u_le_x;
    double velocityY = S->velocity_field[3 * index + 1];
    double velocityZ = S->velocity_field[3 * index + 2];
    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
    // Equation 3.4 with c_s^2 = 1/3
    return weight * S->density_field[index] * (1.0 + dot_product / (S->c_s * S->c_s) + dot_product * dot_product / (2 * S->c_s * S->c_s * S->c_s * S->c_s) - norm_square / (2 * S->c_s * S->c_s));
}

void collision(struct LBM* S) {
    const double tauinv = 1.0 / S->tau;
    const double omtauinv = 1.0 - tauinv;

    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                for(int i = 0; i < S->direction_size; i++) {
                    int feqIndex = (z * S->nX * S->nY) + (y * S->nX) + x;
                    double feq = calculate_feq(feqIndex, S, i);

                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    // Equation 3.9
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tauinv * feq;
                }
            }
        }
    }
}
void periodic_boundary_condition(struct LBM* S) {
    for(int x = 0; x < S->nX; x++) {
        for (int y = 0; y < S->nY; y++) {
            for (int z = 0; z < S->nZ; z++) {
                for (int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    int xmd = (S->nX + x - S->directions[3 * i]) % S->nX;
                    int ymd = (S->nY + y - S->directions[3 * i + 1]) % S->nY;
                    int zmd = (S->nZ + z - S->directions[3 * i + 2]) % S->nZ;
                    int otherIndex = xmd + ymd * S->nX + zmd * S->nX * S->nY + i * S->nX * S->nY * S->nZ;

                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}
void couette_boundary_condition(struct LBM* S) {
    double c_s_square = S->c_s * S->c_s;
    for(int x = 0; x < S->nX; x++) {
        for (int y = 0; y < S->nY; y++) {
            for (int z = 0; z < S->nZ; z++) {
                for (int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    int reverseIndex =  x + y * S->nX + z * S->nX * S->nY + S->reverse_indexes[i] * S->nX * S->nY * S->nZ;
                    int xDirection = S->directions[3 * i];
                    int yDirection = S->directions[3 * i + 1];
                    int zDirection = S->directions[3 * i + 2];
                    if (y == 0 && yDirection == 1) {
                        // Bottom Wall.
                        // Equation 5.27 from LBM Principles and Practice.
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == S->nY - 1 && yDirection == -1) {
                        // Top wall
                        // Equation 5.28 from LBM Principles and Practice.
                        // Coefficients of Equation 5.28 calculated from footnote 17.
                        double u_max = 0.1;
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection * 2 * S->weights[i] / (c_s_square) * u_max;
                    } else {
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

void lees_edwards_boundary_condition(struct LBM* S, int time) {
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
                        int x_pos = (x + d_x_I + 2 * S->nX - directionX) % S->nX;
                        int x_shifted = (x + d_x_I + 2 * S->nX + 1 - directionX) % S->nX;
                        int posIndex = x_pos + ymd * S->nX + zmd * S->nX * S->nY;
                        int shiftedIndex = x_shifted + ymd * S->nX + zmd * S->nX * S->nY;
                        int distIndex = i * S->nX * S->nY * S->nZ;

                        double u_le_x = -1 * S->gamma_dot * (double)S->nY;
                        double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] +
                                calculate_feq_u(posIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(posIndex, 0, S, directionX, directionY, directionZ, S->weights[i]);
                        double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u(shiftedIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(shiftedIndex, 0, S, directionX, directionY, directionZ, S->weights[i]);
                        // Equation (18) from the same paper.
                        S->particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos;
                    } else if(y==S->nY-1 && directionY == -1) {
                        // Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                        // by Eric Lorenz and Alfons G. Hoekstra
                        // Top Wall.
                        int x_pos = (x - d_x_I + 2 * S->nX - directionX) % S->nX;
                        int x_shifted = (x - d_x_I + 2 * S->nX - 1 - directionX) % S->nX;
                        int posIndex = x_pos + ymd * S->nX + zmd * S->nX * S->nY;
                        int shiftedIndex = x_shifted + ymd * S->nX + zmd * S->nX * S->nY;
                        int distIndex = i * S->nX * S->nY * S->nZ;

                        double u_le_x = S->gamma_dot * (double)S->nY;
                        double galilean_transformation_pos = S->previous_particle_distributions[posIndex + distIndex] +
                                calculate_feq_u(posIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(posIndex, 0, S, directionX, directionY, directionZ, S->weights[i]);
                        double galilean_transformation_shift = S->previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u(shiftedIndex, u_le_x, S, directionX, directionY, directionZ, S->weights[i]) -
                                calculate_feq_u(shiftedIndex, 0, S, directionX, directionY, directionZ, S->weights[i]);
                        // Equation (18) from the same paper.
                        S->particle_distributions[index] = s_1 * galilean_transformation_shift + s_2 * galilean_transformation_pos;
                    } else {
                        // Interior points. (Enforce periodicity along x-axis).
                        int xmd = (S->nX + x - directionX) % S->nX;
                        int otherIndex = xmd + ymd * S->nX + zmd * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}
void stream(struct LBM* S, int time) {
    if (S->boundary_condition == 1) {
        periodic_boundary_condition(S);
    } else if (S->boundary_condition == 2) {
        couette_boundary_condition(S);
    } else if (S->boundary_condition == 3) {
        lees_edwards_boundary_condition(S, time);
    }
}

void perform_timestep(struct LBM* S, int time) {
    //fprintf(stderr,"compute_density_momentum_moment \n");
    compute_density_momentum_moment(S);
    //fprintf(stderr,"collision \n");
    collision(S);
    //fprintf(stderr,"stream \n");
    stream(S, time);
}

#endif

// int time, double tau, double gamma_dot, double c_s, int boundary_condition,