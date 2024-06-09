#include "LBM.h"
#include <immintrin.h>


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
                    double dist = S->particle_distributions[index];
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
    {
        int dir_x = directions[0];
        int dir_y = directions[1];
        int dir_z = directions[2];
        for (int index = 0; index < SnXSnYSnZ; index++) {
            double dist = particle_distributions[index];
            density_field[index] = dist;
            velocity_field[3 * index] = dist * dir_x;
            velocity_field[3 * index + 1] = dist * dir_y;
            velocity_field[3 * index + 2] = dist * dir_z;
        }
    }
    for (int i = 1; i < direction_size; i++) {
        int idx = 3 * i;
        int dir_x = directions[idx];
        int dir_y = directions[idx + 1];
        int dir_z = directions[idx + 2];

        for (int index = 0; index < SnXSnYSnZ; index++) {
            double dist = particle_distributions[index + i * SnXSnYSnZ];
            density_field[index] += dist;
            velocity_field[3 * index] += dist * dir_x;
            velocity_field[3 * index + 1] += dist * dir_y;
            velocity_field[3 * index + 2] += dist * dir_z;
        }
    }

    int index;
    for (index = 0; index < SnXSnYSnZ; index++) {
        double new_density = density_field[index];
        double inv_new_density = 1.0 / new_density;

        velocity_field[3 * index] = velocity_field[3 * index] * inv_new_density;
        velocity_field[3 * index + 1] = velocity_field[3 * index + 1] * inv_new_density;
        velocity_field[3 * index + 2] = velocity_field[3 * index + 2] * inv_new_density;
    }
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
    double* velocity_fieldX = S->velocity_fieldX;
    double* velocity_fieldY = S->velocity_fieldY;
    double* velocity_fieldZ = S->velocity_fieldZ;
    {
        int dir_x = directions[0];
        int dir_y = directions[1];
        int dir_z = directions[2];
        for (int index = 0; index < SnXSnYSnZ; index++) {
            double dist = particle_distributions[index];
            density_field[index] = dist;
            velocity_fieldX[index] = dist * dir_x;
            velocity_fieldY[index] = dist * dir_y;
            velocity_fieldZ[index] = dist * dir_z;
        }
    }
    for (int i = 1; i < direction_size; i++) {
        int idx = 3 * i;
        int dir_x = directions[idx];
        int dir_y = directions[idx + 1];
        int dir_z = directions[idx + 2];

        for (int index = 0; index < SnXSnYSnZ; index++) {
            double dist = particle_distributions[index + i * SnXSnYSnZ];
            density_field[index] += dist;
            velocity_fieldX[index] += dist * dir_x;
            velocity_fieldY[index] += dist * dir_y;
            velocity_fieldZ[index] += dist * dir_z;
        }
    }

    // Unrolling the loop by a factor of 4
    int index;
    for (index = 0; index <= SnXSnYSnZ - 4; index += 4) {
        double new_density0 = density_field[index];
        double new_density1 = density_field[index + 1];
        double new_density2 = density_field[index + 2];
        double new_density3 = density_field[index + 3];
        double inv_new_density0 = 1.0 / new_density0;
        double inv_new_density1 = 1.0 / new_density1;
        double inv_new_density2 = 1.0 / new_density2;
        double inv_new_density3 = 1.0 / new_density3;

        velocity_fieldX[index] = velocity_fieldX[index] * inv_new_density0;
        velocity_fieldY[index] = velocity_fieldY[index] * inv_new_density0;
        velocity_fieldZ[index] = velocity_fieldZ[index] * inv_new_density0;

        velocity_fieldX[index + 1] = velocity_fieldX[index + 1] * inv_new_density1;
        velocity_fieldY[index + 1] = velocity_fieldY[index + 1] * inv_new_density1;
        velocity_fieldZ[index + 1] = velocity_fieldZ[index + 1] * inv_new_density1;

        velocity_fieldX[index + 2] = velocity_fieldX[index + 2] * inv_new_density2;
        velocity_fieldY[index + 2] = velocity_fieldY[index + 2] * inv_new_density2;
        velocity_fieldZ[index + 2] = velocity_fieldZ[index + 2] * inv_new_density2;

        velocity_fieldX[index + 3] = velocity_fieldX[index + 3] * inv_new_density3;
        velocity_fieldY[index + 3] = velocity_fieldY[index + 3] * inv_new_density3;
        velocity_fieldZ[index + 3] = velocity_fieldZ[index + 3] * inv_new_density3;
    }

    // Handle the remaining elements
    for (; index < SnXSnYSnZ; index++) {
        double new_density = density_field[index];
        double inv_new_density = 1.0 / new_density;
        velocity_fieldX[index] = velocity_fieldX[index] * inv_new_density;
        velocity_fieldY[index] = velocity_fieldY[index] * inv_new_density;
        velocity_fieldZ[index] = velocity_fieldZ[index] * inv_new_density;
    }
}


void momentumO3(struct LBMarrays* S) {
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


        for (int index = 0; index < SnXSnYSnZ; index += 4) {
            double dist[4];

            // Load dist values for 4 indices in parallel
            for (int j = 0; j < 4; j++) {
                int ind = index + j;
                dist[j] = particle_distributions[ind + i * SnXSnYSnZ];
            }

            // Update temp_density_field for 4 indices in parallel
            __m256d dist_vec = _mm256_loadu_pd(dist);
            __m256d temp_density_field_vec = _mm256_loadu_pd(&temp_density_field[index]);
            temp_density_field_vec = _mm256_add_pd(temp_density_field_vec, dist_vec);
            _mm256_storeu_pd(&temp_density_field[index], temp_density_field_vec);

            // Update temp_velocity_field for 4 indices in parallel
            __m256d dir_x_vec = _mm256_set1_pd(dir_x);
            __m256d dir_y_vec = _mm256_set1_pd(dir_y);
            __m256d dir_z_vec = _mm256_set1_pd(dir_z);
            __m256d dist_dir_x = _mm256_mul_pd(dist_vec, dir_x_vec);
            __m256d dist_dir_y = _mm256_mul_pd(dist_vec, dir_y_vec);
            __m256d dist_dir_z = _mm256_mul_pd(dist_vec, dir_z_vec);

            __m256d temp_velocity_x = _mm256_loadu_pd(&temp_velocity_field[3 * index]);
            __m256d temp_velocity_y = _mm256_loadu_pd(&temp_velocity_field[3 * index + 1]);
            __m256d temp_velocity_z = _mm256_loadu_pd(&temp_velocity_field[3 * index + 2]);

            temp_velocity_x = _mm256_add_pd(temp_velocity_x, dist_dir_x);
            temp_velocity_y = _mm256_add_pd(temp_velocity_y, dist_dir_y);
            temp_velocity_z = _mm256_add_pd(temp_velocity_z, dist_dir_z);

            _mm256_storeu_pd(&temp_velocity_field[3 * index], temp_velocity_x);
            _mm256_storeu_pd(&temp_velocity_field[3 * index + 1], temp_velocity_y);
            _mm256_storeu_pd(&temp_velocity_field[3 * index + 2], temp_velocity_z);
        }
    }

    // Copy the results from temporary arrays to the final arrays
    for (int index = 0; index < SnXSnYSnZ; index += 4) {
        // Load new density values from the temp_density_field
        __m256d new_density_vec = _mm256_loadu_pd(&temp_density_field[index]);

        // Store new density values to density_field
        _mm256_storeu_pd(&density_field[index], new_density_vec);

        // Compute the inverse of new density
        __m256d inv_new_density_vec = _mm256_div_pd(_mm256_set1_pd(1.0), new_density_vec);

        // Load temp velocity values
        __m256d temp_velocity_x_vec = _mm256_loadu_pd(&temp_velocity_field[3 * index]);
        __m256d temp_velocity_y_vec = _mm256_loadu_pd(&temp_velocity_field[3 * index + 4]);
        __m256d temp_velocity_z_vec = _mm256_loadu_pd(&temp_velocity_field[3 * index + 8]);

        // Compute final velocity by multiplying with the inverse density
        __m256d velocity_x_vec = _mm256_mul_pd(temp_velocity_x_vec, inv_new_density_vec);
        __m256d velocity_y_vec = _mm256_mul_pd(temp_velocity_y_vec, inv_new_density_vec);
        __m256d velocity_z_vec = _mm256_mul_pd(temp_velocity_z_vec, inv_new_density_vec);

        // Store the computed velocity values to velocity_field
        _mm256_storeu_pd(&velocity_field[3 * index], velocity_x_vec);
        _mm256_storeu_pd(&velocity_field[3 * index + 4], velocity_y_vec);
        _mm256_storeu_pd(&velocity_field[3 * index + 8], velocity_z_vec);
    }

    // Free temporary arrays
    free(temp_density_field);
    free(temp_velocity_field);
}




void momentumO6(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_fieldX = S->velocity_fieldX;
    double* velocity_fieldY = S->velocity_fieldY;
    double* velocity_fieldZ = S->velocity_fieldZ;

    {
        __m256d vDirectionX = _mm256_set1_pd(directions[0]);
        __m256d vDirectionY = _mm256_set1_pd(directions[1]);
        __m256d vDirectionZ = _mm256_set1_pd(directions[2]);
        for (int index = 0; index < SnXSnYSnZ; index+=4) {
            __m256d distV = _mm256_loadu_pd(&particle_distributions[index]);
            _mm256_storeu_pd(&density_field[index], distV);
            _mm256_storeu_pd(&velocity_fieldX[index], _mm256_mul_pd(distV, vDirectionX));
            _mm256_storeu_pd(&velocity_fieldY[index], _mm256_mul_pd(distV, vDirectionY));
            _mm256_storeu_pd(&velocity_fieldZ[index], _mm256_mul_pd(distV, vDirectionZ));
        }
    }
    for (int i = 1; i < direction_size; i++) {
        int idx = 3 * i;
        __m256d vDirectionX = _mm256_set1_pd(directions[idx]);
        __m256d vDirectionY = _mm256_set1_pd(directions[idx + 1]);
        __m256d vDirectionZ = _mm256_set1_pd(directions[idx + 2]);
        int particleIndex = i * SnXSnYSnZ;
        for (int index = 0; index < SnXSnYSnZ; index+=4) {
            __m256d distV = _mm256_loadu_pd(&particle_distributions[index + particleIndex]);
            _mm256_storeu_pd(&density_field[index], _mm256_add_pd(_mm256_loadu_pd(&density_field[index]), distV));
            _mm256_storeu_pd(&velocity_fieldX[index], _mm256_fmadd_pd(distV, vDirectionX, _mm256_loadu_pd(&velocity_fieldX[index])));
            _mm256_storeu_pd(&velocity_fieldY[index], _mm256_fmadd_pd(distV, vDirectionY, _mm256_loadu_pd(&velocity_fieldY[index])));
            _mm256_storeu_pd(&velocity_fieldZ[index], _mm256_fmadd_pd(distV, vDirectionZ, _mm256_loadu_pd(&velocity_fieldZ[index])));
        }
    }

    // Unrolling the loop by a factor of 4
    __m256d one = _mm256_set1_pd(1.0);
    for (int index = 0; index < SnXSnYSnZ; index += 4) {
        __m256d density = _mm256_loadu_pd(&density_field[index]);
        __m256d inv_density = _mm256_div_pd(one, density);
        __m256d vX = _mm256_loadu_pd(&velocity_fieldX[index]);
        _mm256_storeu_pd(&velocity_fieldX[index], _mm256_mul_pd(vX, inv_density));
        __m256d vY = _mm256_loadu_pd(&velocity_fieldY[index]);
        _mm256_storeu_pd(&velocity_fieldY[index], _mm256_mul_pd(vY, inv_density));
        __m256d vZ = _mm256_loadu_pd(&velocity_fieldZ[index]);
        _mm256_storeu_pd(&velocity_fieldZ[index], _mm256_mul_pd(vZ, inv_density));
    }
}


void momentumO61(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* density_field = S->density_field;
    double* velocity_fieldX = S->velocity_fieldX;
    double* velocity_fieldY = S->velocity_fieldY;
    double* velocity_fieldZ = S->velocity_fieldZ;

    {
        __m256d vDirectionX = _mm256_set1_pd(directions[0]);
        __m256d vDirectionY = _mm256_set1_pd(directions[1]);
        __m256d vDirectionZ = _mm256_set1_pd(directions[2]);
        for (int index = 0; index < SnXSnYSnZ; index+=4) {
            __m256d distV = _mm256_loadu_pd(&particle_distributions[index]);
            _mm256_storeu_pd(&density_field[index], distV);
            __m256d AX =_mm256_mul_pd(distV, vDirectionX);
            __m256d AY =_mm256_mul_pd(distV, vDirectionY);
            __m256d AZ =_mm256_mul_pd(distV, vDirectionZ);
            _mm256_storeu_pd(&velocity_fieldX[index], AX);
            _mm256_storeu_pd(&velocity_fieldY[index], AY);
            _mm256_storeu_pd(&velocity_fieldZ[index], AZ);
        }
    }
    for (int i = 1; i < direction_size; i++) {
        int idx = 3 * i;
        __m256d vDirectionX = _mm256_set1_pd(directions[idx]);
        __m256d vDirectionY = _mm256_set1_pd(directions[idx + 1]);
        __m256d vDirectionZ = _mm256_set1_pd(directions[idx + 2]);
        int particleIndex = i * SnXSnYSnZ;
        for (int index = 0; index < SnXSnYSnZ; index+=4) {
            __m256d distV = _mm256_loadu_pd(&particle_distributions[index + particleIndex]);
            __m256d A = _mm256_loadu_pd(&density_field[index]);
            __m256d B = _mm256_add_pd(A, distV);
            _mm256_storeu_pd(&density_field[index], B);
            __m256d AX = _mm256_loadu_pd(&velocity_fieldX[index]);
            __m256d AY = _mm256_loadu_pd(&velocity_fieldY[index]);
            __m256d AZ = _mm256_loadu_pd(&velocity_fieldZ[index]);     
            __m256d BX = _mm256_fmadd_pd(distV, vDirectionX, AX);            
            __m256d BY = _mm256_fmadd_pd(distV, vDirectionY, AY);
            __m256d BZ = _mm256_fmadd_pd(distV, vDirectionZ, AZ);
            _mm256_storeu_pd(&velocity_fieldX[index], BX);
            _mm256_storeu_pd(&velocity_fieldY[index], BY);
            _mm256_storeu_pd(&velocity_fieldZ[index], BZ);
        }
    }

    // Unrolling the loop by a factor of 4
    __m256d one = _mm256_set1_pd(1.0);
    for (int index = 0; index < SnXSnYSnZ; index += 4) {
        __m256d density = _mm256_loadu_pd(&density_field[index]);
        __m256d inv_density = _mm256_div_pd(one, density);
        __m256d vX = _mm256_loadu_pd(&velocity_fieldX[index]);
        __m256d vY = _mm256_loadu_pd(&velocity_fieldY[index]);
        __m256d vZ = _mm256_loadu_pd(&velocity_fieldZ[index]);
        __m256d AX =  _mm256_mul_pd(vX, inv_density);
        __m256d AY =  _mm256_mul_pd(vY, inv_density);
        __m256d AZ =  _mm256_mul_pd(vZ, inv_density);
        _mm256_storeu_pd(&velocity_fieldX[index], AX);
        _mm256_storeu_pd(&velocity_fieldY[index], AY);
        _mm256_storeu_pd(&velocity_fieldZ[index], AZ);
    }
}



////--- OPTIMIZATION STEP ONE - AVX ----- //////////////

/**
 * These functions calculate the macroscopic moments and the equilibrium distributions as described in
 * the steps 1 & 2 of the Time Step Algorithm described in Section 3.3.2 of The Lattice Boltzmann Method book.
 */

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