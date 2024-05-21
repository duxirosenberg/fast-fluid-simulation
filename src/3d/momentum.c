#include "LBM.h"
//#include <immintrin.h>
#include <stdlib.h>


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



////-------- OPTIMIZATION STEP ONE ------//////////////



// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 3 + nX * nY (1 + nZ * (3 + 4 * q))

void momentumO1(struct LBMarrays* S) {

    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    int ySnX, zSnXSnY, index, idx;
    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            ySnX = y * SnX; 
            for(int z = 0; z < S->nZ; z++) {
                zSnXSnY = z * SnXSnY;
                double new_density = 0;
                double u[] = {0, 0, 0};
                for(int i = 0; i < S->direction_size; i++) {


                    index = x + ySnX + zSnXSnY + i * SnXSnYSnZ;
                    double dist = S->particle_distributions[index];;
                    new_density += dist;
                    idx = 3 * i;
                    u[0] += dist * (double)S->directions[idx];
                    u[1] += dist * (double)S->directions[idx + 1];
                    u[2] += dist * (double)S->directions[idx + 2];
                }
                index = zSnXSnY + ySnX + x;
                S->density_field[index] = new_density;
                double inv_new_density = 1.0 / new_density;
                velocity_field[3 * index] = u[0] * inv_new_density;
                velocity_field[3 * index + 1] = u[1] * inv_new_density;
                velocity_field[3 * index + 2] = u[2] * inv_new_density;
            }
        }
    }
}


////-------- OPTIMIZATION STEP TWO -------- //////////////


// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 2 + q + 4 * q * nX * nY * nZ
void momentumO2(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    // Temporary arrays for accumulation
    double* temp_density_field = (double*)calloc(SnXSnYSnZ, sizeof(double));
    double* temp_velocity_field = (double*)calloc(3 * SnXSnYSnZ, sizeof(double));

    for (int i = 0; i < direction_size; i++) {
        int idx = 3 * i;
        int dir_x = directions[idx];
        int dir_y = directions[idx + 1];
        int dir_z = directions[idx + 2];

        for (int index = 0; index < SnXSnYSnZ; index++) {
            int x = index % SnX;
            int y = (index / SnX) % SnY;
            int z = index / SnXSnY;

            double dist = particle_distributions[index + i * SnXSnYSnZ];
            temp_density_field[index] += dist;
            temp_velocity_field[3 * index] += dist * dir_x;
            temp_velocity_field[3 * index + 1] += dist * dir_y;
            temp_velocity_field[3 * index + 2] += dist * dir_z;
        }
    }

    // Copy the results from temporary arrays to the final arrays
    for (int index = 0; index < SnXSnYSnZ; index++) {
        double new_density = temp_density_field[index];
        density_field[index] = new_density;

        double inv_new_density = 1.0 / new_density;
        velocity_field[3 * index] = temp_velocity_field[3 * index] * inv_new_density;
        velocity_field[3 * index + 1] = temp_velocity_field[3 * index + 1] * inv_new_density;
        velocity_field[3 * index + 2] = temp_velocity_field[3 * index + 2] * inv_new_density;
    }

    // Free temporary arrays
    free(temp_density_field);
    free(temp_velocity_field);
}


// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 2 + q + 4 * q * nX * nY * nZ

void momentumO25(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    // Temporary arrays for accumulation
    double* temp_density_field = (double*)calloc(SnXSnYSnZ, sizeof(double));
    double* temp_velocity_field = (double*)calloc(3 * SnXSnYSnZ, sizeof(double));

    for (int i = 0; i < direction_size; i++) {
        int idx = 3 * i;
        int dir_x = directions[idx];
        int dir_y = directions[idx + 1];
        int dir_z = directions[idx + 2];

        for (int index = 0; index < SnXSnYSnZ; index++) {
            int x = index % SnX;
            int y = (index / SnX) % SnY;
            int z = index / SnXSnY;

            double dist = particle_distributions[index + i * SnXSnYSnZ];
            temp_density_field[index] += dist;
            temp_velocity_field[3 * index] += dist * dir_x;
            temp_velocity_field[3 * index + 1] += dist * dir_y;
            temp_velocity_field[3 * index + 2] += dist * dir_z;
        }
    }

    // Unrolling the loop by a factor of 4
    int index;
    for (index = 0; index <= SnXSnYSnZ - 4; index += 4) {
        double new_density0 = temp_density_field[index];
        double new_density1 = temp_density_field[index + 1];
        double new_density2 = temp_density_field[index + 2];
        double new_density3 = temp_density_field[index + 3];

        density_field[index] = new_density0;
        density_field[index + 1] = new_density1;
        density_field[index + 2] = new_density2;
        density_field[index + 3] = new_density3;

        double inv_new_density0 = 1.0 / new_density0;
        double inv_new_density1 = 1.0 / new_density1;
        double inv_new_density2 = 1.0 / new_density2;
        double inv_new_density3 = 1.0 / new_density3;

        velocity_field[3 * index] = temp_velocity_field[3 * index] * inv_new_density0;
        velocity_field[3 * index + 1] = temp_velocity_field[3 * index + 1] * inv_new_density0;
        velocity_field[3 * index + 2] = temp_velocity_field[3 * index + 2] * inv_new_density0;

        velocity_field[3 * (index + 1)] = temp_velocity_field[3 * (index + 1)] * inv_new_density1;
        velocity_field[3 * (index + 1) + 1] = temp_velocity_field[3 * (index + 1) + 1] * inv_new_density1;
        velocity_field[3 * (index + 1) + 2] = temp_velocity_field[3 * (index + 1) + 2] * inv_new_density1;

        velocity_field[3 * (index + 2)] = temp_velocity_field[3 * (index + 2)] * inv_new_density2;
        velocity_field[3 * (index + 2) + 1] = temp_velocity_field[3 * (index + 2) + 1] * inv_new_density2;
        velocity_field[3 * (index + 2) + 2] = temp_velocity_field[3 * (index + 2) + 2] * inv_new_density2;

        velocity_field[3 * (index + 3)] = temp_velocity_field[3 * (index + 3)] * inv_new_density3;
        velocity_field[3 * (index + 3) + 1] = temp_velocity_field[3 * (index + 3) + 1] * inv_new_density3;
        velocity_field[3 * (index + 3) + 2] = temp_velocity_field[3 * (index + 3) + 2] * inv_new_density3;
    }

    // Handle the remaining elements
    for (; index < SnXSnYSnZ; index++) {
        double new_density = temp_density_field[index];
        density_field[index] = new_density;

        double inv_new_density = 1.0 / new_density;
        velocity_field[3 * index] = temp_velocity_field[3 * index] * inv_new_density;
        velocity_field[3 * index + 1] = temp_velocity_field[3 * index + 1] * inv_new_density;
        velocity_field[3 * index + 2] = temp_velocity_field[3 * index + 2] * inv_new_density;
    }

    // Free temporary arrays
    free(temp_density_field);
    free(temp_velocity_field);
}



////--- OPTIMIZATION STEP ONE - AVX ----- //////////////

/**
 * These functions calculate the macroscopic moments and the equilibrium distributions as described in
 * the steps 1 & 2 of the Time Step Algorithm described in Section 3.3.2 of The Lattice Boltzmann Method book.
 */


/*void momentumO6(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    for (int z = 0; z < SnZ; z++) {
        int zSnXSnY = z * SnXSnY;
        for (int y = 0; y < SnY; y++) {
            int ySnX = y * SnX;
            for (int x = 0; x < SnX; x += 4) { // Process 4 elements at a time using AVX256
                int index = x + ySnX + zSnXSnY;
                double new_density[4] = {0}; // Store 4 densities
                double u[3][4] = {{0}}; // Store 3 velocities for 4 elements

                for (int i = 0; i < direction_size; i++) {
                    double dist[4];
                    int idx = 3 * i;
                    __m256d direction_x = _mm256_set1_pd(directions[idx]);
                    __m256d direction_y = _mm256_set1_pd(directions[idx + 1]);
                    __m256d direction_z = _mm256_set1_pd(directions[idx + 2]);

                    // Load 4 particle_distributions elements
                    __m256d pd = _mm256_loadu_pd(&particle_distributions[index]);

                    // Calculate 4 new_densities
                    __m256d nd = _mm256_add_pd(_mm256_loadu_pd(new_density), pd);
                    _mm256_storeu_pd(new_density, nd);

                    // Calculate 4 u[0] components
                    __m256d ud_x = _mm256_mul_pd(pd, direction_x);
                    _mm256_storeu_pd(u[0], _mm256_add_pd(_mm256_loadu_pd(u[0]), ud_x));

                    // Calculate 4 u[1] components
                    __m256d ud_y = _mm256_mul_pd(pd, direction_y);
                    _mm256_storeu_pd(u[1], _mm256_add_pd(_mm256_loadu_pd(u[1]), ud_y));

                    // Calculate 4 u[2] components
                    __m256d ud_z = _mm256_mul_pd(pd, direction_z);
                    _mm256_storeu_pd(u[2], _mm256_add_pd(_mm256_loadu_pd(u[2]), ud_z));

                    index += SnXSnYSnZ;
                }

                index = x + ySnX + zSnXSnY;
                for (int i = 0; i < 4; i++) {
                    density_field[index + i] = new_density[i];
                    double inv_new_density = 1.0 / new_density[i];
                    velocity_field[3 * (index + i)] = u[0][i] * inv_new_density;
                    velocity_field[3 * (index + i) + 1] = u[1][i] * inv_new_density;
                    velocity_field[3 * (index + i) + 2] = u[2][i] * inv_new_density;
                }
            }
        }
    }
}*/


/*// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: 3 + nX * nY (1 + nZ * (5 + q))
void momentumO11 (struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnX * SnY * SnZ;

    int ySnX, zSnXSnY, index;


    for(int x = 0; x < SnX; x++) {
        for(int y = 0; y < SnY; y++) {

            ySnX = y * SnX; 
            for(int z = 0; z < SnZ; z++) {

                zSnXSnY = z * SnXSnY;
                double new_density = 0;
                double u[] = {0, 0, 0};
                index = x + ySnX + zSnXSnY;
                for(int i = 0; i < S->direction_size; i++) {
                    
                    double dist = S->particle_distributions[index];
                    new_density += dist;
                    u[0] += dist * (double)S->directions[3 * i];
                    u[1] += dist * (double)S->directions[3 * i + 1];
                    u[2] += dist * (double)S->directions[3 * i + 2];

                    index += SnXSnYSnZ;
                }
                index = zSnXSnY + ySnX + x;
                S->density_field[index] = new_density;
                S->velocity_field[3 * index] = u[0] / new_density;
                S->velocity_field[3 * index + 1] = u[1] / new_density;
                S->velocity_field[3 * index + 2] = u[2] / new_density;
            }
        }
    }
}

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: 3 + nX * nY (1 + nZ * (3 + 4 * q))
void momentumO12(struct LBMarrays* S) {

    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnX * SnY * SnZ;

    int ySnX, zSnXSnY, index;
    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            ySnX = y * SnX; 
            for(int z = 0; z < S->nZ; z++) {
                zSnXSnY = z * SnXSnY;
                double new_density = 0;
                double u[] = {0, 0, 0};
                for(int i = 0; i < S->direction_size; i++) {


                    int index = x + ySnX + zSnXSnY + i * SnXSnYSnZ;
                    double dist = S->particle_distributions[index];;
                    new_density += dist;
                    u[0] += dist * (double)S->directions[3 * i];
                    u[1] += dist * (double)S->directions[3 * i + 1];
                    u[2] += dist * (double)S->directions[3 * i + 2];
                }
                index = zSnXSnY + ySnX + x;
                S->density_field[index] = new_density;
                S->velocity_field[3 * index] = u[0] / new_density;
                S->velocity_field[3 * index + 1] = u[1] / new_density;
                S->velocity_field[3 * index + 2] = u[2] / new_density;
            }
        }
    }
}*/


/*//loop order Z - Y - X

// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 3 + nZ(1 + nY(1 + nX(4 + 5q))) 
void momentumO24(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    for (int z = 0; z < SnZ; z++) {
        int zSnXSnY = z * SnXSnY;
        for (int y = 0; y < SnY; y++) {
            int ySnX = y * SnX;
            for (int x = 0; x < SnX; x++) {
                int index = x + ySnX + zSnXSnY;
                double new_density = 0;
                double u[3] = {0, 0, 0};

                for (int i = 0; i < direction_size; i++) {
                    double dist = particle_distributions[index];
                    new_density += dist;
                    int idx = 3 * i;
                    u[0] += dist * directions[idx];
                    u[1] += dist * directions[idx + 1];
                    u[2] += dist * directions[idx + 2];
                    index += SnXSnYSnZ;
                }

                index = x + ySnX + zSnXSnY;
                density_field[index] = new_density;
                double inv_new_density = 1.0 / new_density;
                velocity_field[3 * index] = u[0] * inv_new_density;
                velocity_field[3 * index + 1] = u[1] * inv_new_density;
                velocity_field[3 * index + 2] = u[2] * inv_new_density;
            }
        }
    }
}





void momentum_arraysO2(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions) {

    int nXnY = nX * nY;
    int nXnYnZ = nXnY * nZ;
    for(int z = 0; z < nZ; z++) {
        int znXnY = z * nXnY;
        for(int y = 0; y < nY; y++) {
            int ySnX = y * nX;
            for(int x = 0; x < nX; x++) {
                //Equation 3.1
                int index = x + ySnX + znXnY;
                double new_density = 0;
                double u[] = {0, 0, 0};
                for(int i = 0; i < direction_size; i++) {
                    double dist = particle_distributions[index];
                    new_density += dist;
                    int idx = 3 * i;
                    u[0] += dist * directions[idx];
                    u[1] += dist * directions[idx + 1];
                    u[2] += dist * directions[idx + 2];
                    index += nXnYnZ;
                }
                index = x + ySnX + znXnY;
                density_field[index] = new_density;
                double inv_new_density = 1.0/ new_density;
                velocity_field[3 * index] = u[0] * inv_new_density;
                velocity_field[3 * index + 1] = u[1] * inv_new_density;
                velocity_field[3 * index + 2] = u[2] * inv_new_density;
            }
        }
    }
}


void momentumO21(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    for (int z = 0; z < SnZ; z++) {
        int zSnXSnY = z * SnXSnY;
        for (int y = 0; y < SnY; y++) {
            int ySnX = y * SnX;
            for (int x = 0; x < SnX; x++) {
                int index = x + ySnX + zSnXSnY;
                double new_density = 0;
                double u[3] = {0, 0, 0};

                int i;
                for (i = 0; i <= direction_size - 4; i += 4) {
                    int idx0 = 3 * i;
                    int idx1 = 3 * (i + 1);
                    int idx2 = 3 * (i + 2);
                    int idx3 = 3 * (i + 3);

                    double dist0 = particle_distributions[index];
                    double dist1 = particle_distributions[index + SnXSnYSnZ];
                    double dist2 = particle_distributions[index + 2 * SnXSnYSnZ];
                    double dist3 = particle_distributions[index + 3 * SnXSnYSnZ];

                    new_density += dist0 + dist1 + dist2 + dist3;

                    u[0] += dist0 * directions[idx0] + dist1 * directions[idx1] + dist2 * directions[idx2] + dist3 * directions[idx3];
                    u[1] += dist0 * directions[idx0 + 1] + dist1 * directions[idx1 + 1] + dist2 * directions[idx2 + 1] + dist3 * directions[idx3 + 1];
                    u[2] += dist0 * directions[idx0 + 2] + dist1 * directions[idx1 + 2] + dist2 * directions[idx2 + 2] + dist3 * directions[idx3 + 2];

                    index += 4 * SnXSnYSnZ;
                }

                // Handle any remaining elements
                for (; i < direction_size; i++) {
                    double dist = particle_distributions[index];
                    new_density += dist;
                    int idx = 3 * i;
                    u[0] += dist * directions[idx];
                    u[1] += dist * directions[idx + 1];
                    u[2] += dist * directions[idx + 2];
                    index += SnXSnYSnZ;
                }

                index = x + y * SnX + z * SnXSnY;
                density_field[index] = new_density;
                double inv_new_density = 1.0 / new_density;
                velocity_field[3 * index] = u[0] * inv_new_density;
                velocity_field[3 * index + 1] = u[1] * inv_new_density;
                velocity_field[3 * index + 2] = u[2] * inv_new_density;
            }
        }
    }
}


void momentumO22(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;


    for (int z = 0; z < SnZ; z++) {
        for (int y = 0; y < SnY; y++) {
            for (int x = 0; x < SnX; x++) {
                int index = x + y * SnX + z * SnXSnY;
                double new_density = 0;
                double u[3] = {0, 0, 0};

                int i;
                for (i = 0; i <= direction_size - 8; i += 8) {
                    int idx0 = 3 * i;
                    int idx1 = 3 * (i + 1);
                    int idx2 = 3 * (i + 2);
                    int idx3 = 3 * (i + 3);
                    int idx4 = 3 * (i + 4);
                    int idx5 = 3 * (i + 5);
                    int idx6 = 3 * (i + 6);
                    int idx7 = 3 * (i + 7);

                    double dist0 = particle_distributions[index];
                    double dist1 = particle_distributions[index + SnXSnYSnZ];
                    double dist2 = particle_distributions[index + 2 * SnXSnYSnZ];
                    double dist3 = particle_distributions[index + 3 * SnXSnYSnZ];
                    double dist4 = particle_distributions[index + 4 * SnXSnYSnZ];
                    double dist5 = particle_distributions[index + 5 * SnXSnYSnZ];
                    double dist6 = particle_distributions[index + 6 * SnXSnYSnZ];
                    double dist7 = particle_distributions[index + 7 * SnXSnYSnZ];

                    new_density += dist0 + dist1 + dist2 + dist3 + dist4 + dist5 + dist6 + dist7;

                    u[0] += dist0 * directions[idx0] + dist1 * directions[idx1] + dist2 * directions[idx2] + dist3 * directions[idx3]
                          + dist4 * directions[idx4] + dist5 * directions[idx5] + dist6 * directions[idx6] + dist7 * directions[idx7];
                    u[1] += dist0 * directions[idx0 + 1] + dist1 * directions[idx1 + 1] + dist2 * directions[idx2 + 1] + dist3 * directions[idx3 + 1]
                          + dist4 * directions[idx4 + 1] + dist5 * directions[idx5 + 1] + dist6 * directions[idx6 + 1] + dist7 * directions[idx7 + 1];
                    u[2] += dist0 * directions[idx0 + 2] + dist1 * directions[idx1 + 2] + dist2 * directions[idx2 + 2] + dist3 * directions[idx3 + 2]
                          + dist4 * directions[idx4 + 2] + dist5 * directions[idx5 + 2] + dist6 * directions[idx6 + 2] + dist7 * directions[idx7 + 2];

                    index += 8 * SnXSnYSnZ;
                }

                // Handle any remaining elements
                for (; i < direction_size; i++) {
                    double dist = particle_distributions[index];
                    new_density += dist;
                    int idx = 3 * i;
                    u[0] += dist * directions[idx];
                    u[1] += dist * directions[idx + 1];
                    u[2] += dist * directions[idx + 2];
                    index += SnXSnYSnZ;
                }

                index = x + y * SnX + z * SnXSnY;
                density_field[index] = new_density;
                double inv_new_density = 1.0 / new_density;
                velocity_field[3 * index] = u[0] * inv_new_density;
                velocity_field[3 * index + 1] = u[1] * inv_new_density;
                velocity_field[3 * index + 2] = u[2] * inv_new_density;
            }
        }
    }
}


// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 3 + nZ(1 + nY(1 + nX(4 + 5q))) 
void momentumO23(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_field = S->velocity_field;

    for (int z = 0; z < SnZ; z++) {
        for (int y = 0; y < SnY; y++) {
            for (int x = 0; x < SnX; x++) {
                int index = x + y * SnX + z * SnXSnY;
                double new_density = 0.0;
                double u[3] = {0.0, 0.0, 0.0};

                for (int i = 0; i < direction_size; i++) {
                    double dist = particle_distributions[index + i * SnXSnYSnZ];
                    new_density += dist;
                    int idx = 3 * i;
                    u[0] += dist * directions[idx];
                    u[1] += dist * directions[idx + 1];
                    u[2] += dist * directions[idx + 2];
                }

                density_field[index] = new_density;
                double inv_new_density = 1.0 / new_density;
                velocity_field[3 * index] = u[0] * inv_new_density;
                velocity_field[3 * index + 1] = u[1] * inv_new_density;
                velocity_field[3 * index + 2] = u[2] * inv_new_density;
            }
        }
    }
}*/