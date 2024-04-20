#include "LBM.h"

#ifndef LBM_STRUCT

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

void compute_density_momentum_moment(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions) {
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

double calculate_feq(int index,
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

double calculate_feq_u(int index,
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

void collision(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
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
                    double feq = calculate_feq(feqIndex, c_s, density_field, velocity_field, (double) directions[3 * i], (double) directions[3 * i + 1], (double) directions[3 * i + 2], weights[i]);

                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    //Equation 3.9
                    previous_particle_distributions[index] = omtauinv * particle_distributions[index] + tauinv * feq;

                }
            }
        }
    }
}

void periodic_boundary_condition(int nX, int nY, int nZ, int direction_size,
                               double* previous_particle_distributions,
                               double* particle_distributions,
                               const int* directions
                               ) {
    for(int x = 0; x < nX; x++) {
        for (int y = 0; y < nY; y++) {
            for (int z = 0; z < nZ; z++) {
                for (int i = 0; i < direction_size; i++) {
                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    int xmd = (nX + x - directions[3 * i]) % nX;
                    int ymd = (nY + y - directions[3 * i + 1]) % nY;
                    int zmd = (nZ + z - directions[3 * i + 2]) % nZ;
                    int otherIndex = xmd + ymd * nX + zmd * nX * nY + i * nX * nY * nZ;

                    particle_distributions[index] = previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}

void couette_boundary_condition(int nX, int nY, int nZ, int direction_size, double c_s,
                              double* previous_particle_distributions,
                              double* particle_distributions,
                              const int* directions,
                              const double* weights,
                              int* reverse_indexes) {
    double c_s_square = c_s * c_s;
    for(int x = 0; x < nX; x++) {
        for (int y = 0; y < nY; y++) {
            for (int z = 0; z < nZ; z++) {
                for (int i = 0; i < direction_size; i++) {
                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    int reverseIndex =  x + y * nX + z * nX * nY + reverse_indexes[i] * nX * nY * nZ;
                    int xDirection = directions[3 * i];
                    int yDirection = directions[3 * i + 1];
                    int zDirection = directions[3 * i + 2];
                    if (y==0 && yDirection == 1) {
                        //Bottom Wall.
                        //Equation 5.27 from LBM Principles and Practice.
                        particle_distributions[index]=previous_particle_distributions[reverseIndex];
                    } else if (y==nY - 1 && yDirection == -1) {
                        //Top wall
                        //Equation 5.28 from LBM Principles and Practice.
                        //coefficients of Equation 5.28 calculated from footnote 17.
                        double u_max = 0.1;
                        particle_distributions[index] = previous_particle_distributions[reverseIndex] + xDirection * 2 * weights[i] / (c_s_square) * u_max;
                        //particle_distributions[scalar_index(x,y,z,(4-1))]=previous_particle_distributions[scalar_index(x,y,z,(2-1))];
                        //particle_distributions[scalar_index(x,y,z,(7-1))]=previous_particle_distributions[scalar_index(x,y,z,(5-1))]-(1.0/6.0)*u_max;
                        //particle_distributions[scalar_index(x,y,z,(8-1))]=previous_particle_distributions[scalar_index(x,y,z,(6-1))]+(1.0/6.0)*u_max;
                    } else {
                        //Chapter 13 Taylor-Green periodicity from same book.
                        int xmd = (nX + x - xDirection) % nX;

                        //int ymd = (NY + y - (int)directions[i].y) % NY;
                        int ymd = y - yDirection;
                        int zmd = (nZ + z - zDirection) % nZ;
                        int otherIndex = xmd + ymd * nX + zmd * nX * nY + i * nX * nY * nZ;
                        particle_distributions[index] = previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}

void lees_edwards_boundary_condition(int nX, int nY, int nZ, int direction_size, int time, double gamma_dot, double c_s,
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
                                calculate_feq_u(posIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u(posIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
                        double galilean_transformation_shift = previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u(shiftedIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u(shiftedIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
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
                                calculate_feq_u(posIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u(posIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
                        double galilean_transformation_shift = previous_particle_distributions[shiftedIndex + distIndex] +
                                calculate_feq_u(shiftedIndex, u_le_x, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]) -
                                calculate_feq_u(shiftedIndex, 0, c_s, density_field, velocity_field, directionX, directionY, directionZ, weights[i]);
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

void stream(int nX, int nY, int nZ, int direction_size, int time, double gamma_dot, double c_s, int boundary_condition,
            double* density_field,
            double* velocity_field,
            double* previous_particle_distributions,
            double* particle_distributions,
            const int* directions,
            const double* weights,
            int* reverse_indexes
            ) {
    if (boundary_condition == 1) {
        periodic_boundary_condition(nX, nY, nZ, direction_size, previous_particle_distributions, particle_distributions, directions);
    } else if(boundary_condition == 2) {
        couette_boundary_condition(nX, nY, nZ, direction_size, c_s, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);
    } else if(boundary_condition == 3) {
        lees_edwards_boundary_condition(nX, nY, nZ, direction_size, time, gamma_dot, c_s, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights);
    }
}

void perform_timestep(int nX, int nY, int nZ, int direction_size, int time, double tau, double gamma_dot, double c_s, int boundary_condition,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights,
                      int* reverse_indexes
) {
    //fprintf(stderr,"collision \n");
    collision(nX, nY, nZ, direction_size, tau, c_s, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights);

    //fprintf(stderr,"stream \n");
    stream(nX, nY, nZ, direction_size, time, gamma_dot, c_s, boundary_condition, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);

    //fprintf(stderr,"compute_density_momentum_moment \n");
    compute_density_momentum_moment(nX, nY, nZ, direction_size, density_field, velocity_field, particle_distributions, directions);
}

#endif