#include "LBM.h"

#define CACHE_BLOCK_SIZE 128  // Define the block size in B for cache optimization


// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
void stream_couette_baseline(struct LBMarrays* S) {
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


/*
  Optimization 1: Code Motion
  ---------------------------
  - Pre-compute and move constant expressions out of loops to reduce redundant calculations.
  - Reuse calculated values wherever possible.
  - Simplify and reorganize the code for clarity and efficiency.
*/
void stream_couette_opt1(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;
    double inv_c_s_square = 1.0 / c_s_square;  // Pre-compute inverse of c_s_square
    double u_max = 0.1;  // Maximum velocity
    
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    
    int nXY = nX * nY;
    int nXYZ = nXY * nZ;
    
    for (int x = 0; x < nX; x++) {
        for (int y = 0; y < nY; y++) {
            for (int z = 0; z < nZ; z++) {
                int baseIndex = x + y * nX + z * nXY; // Precompute base index
                
                for (int i = 0; i < q; i++) {
                    int index = baseIndex + i * nXYZ; // Use precomputed base index
                    int reverseIndex = baseIndex + S->reverse_indexes[i] * nXYZ;
                    int xDirection = S->directions[3 * i];
                    int yDirection = S->directions[3 * i + 1];
                    int zDirection = S->directions[3 * i + 2];
                    
                    if (y == 0 && yDirection == 1) {
                        // Bottom Wall.
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == nY - 1 && yDirection == -1) {
                        // Top wall
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection * 2 * S->weights[i] * inv_c_s_square * u_max;
                    } else {
                        // Taylor-Green periodicity
                        int xmd = (nX + x - xDirection) % nX;
                        int ymd = y - yDirection;
                        int zmd = (nZ + z - zDirection) % nZ;
                        int otherIndex = xmd + ymd * nX + zmd * nXY + i * nXYZ;
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}



/*
  Optimization 2: Loop Structure and Blocking
  -------------------------------------------
  - Introduce cache blocking to improve cache locality and reduce cache misses.
  - Reorganize loops to iterate over blocks of data that fit within the cache.
  - Maintain the benefits of code motion from Optimization 1.
*/
void stream_couette_opt2(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;
    double inv_c_s_square = 1.0 / c_s_square;  // Pre-compute inverse of c_s_square
    double u_max = 0.1;  // Maximum velocity

    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;

    int nXY = S->nXY;
    int nXYZ = S->nXYZ;
    int blockSize = CACHE_BLOCK_SIZE / sizeof(double);  // Convert cache block size to the number of doubles

    // Loop over all directions first for better cache locality
    for (int i = 0; i < q; i++) {
        int xDirection = S->directions[3 * i];
        int yDirection = S->directions[3 * i + 1];
        int zDirection = S->directions[3 * i + 2];

        for (int xb = 0; xb < nX; xb += blockSize) {
            for (int yb = 0; yb < nY; yb += blockSize) {
                for (int zb = 0; zb < nZ; zb += blockSize) {
                    int xMax = (xb + blockSize < nX) ? xb + blockSize : nX;
                    int yMax = (yb + blockSize < nY) ? yb + blockSize : nY;
                    int zMax = (zb + blockSize < nZ) ? zb + blockSize : nZ;

                    for (int z = zb; z < zMax; z++) {
                        for (int y = yb; y < yMax; y++) {
                            for (int x = xb; x < xMax; x++) {
                                int baseIndex = x + y * nX + z * nXY; // Precompute base index
                                int index = baseIndex + i * nXYZ; // Use precomputed base index
                                int reverseIndex = baseIndex + S->reverse_indexes[i] * nXYZ;

                                if (y == 0 && yDirection == 1) {
                                    // Bottom Wall.
                                    S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                                } else if (y == nY - 1 && yDirection == -1) {
                                    // Top wall
                                    S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection * 2 * S->weights[i] * inv_c_s_square * u_max;
                                } else {
                                    // Taylor-Green periodicity
                                    int xmd = (nX + x - xDirection) % nX;
                                    int ymd = y - yDirection;
                                    int zmd = (nZ + z - zDirection) % nZ;
                                    int otherIndex = xmd + ymd * nX + zmd * nXY + i * nXYZ;
                                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}








// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
void stream_couette_arrays(int nX, int nY, int nZ, int direction_size, double c_s,
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
                        particle_distributions[index]=previous_particle_distributions[reverseIndex];
                    } else if (y==nY - 1 && yDirection == -1) {
                        double u_max = 0.1;
                        particle_distributions[index] = previous_particle_distributions[reverseIndex] + xDirection * 2 * weights[i] / (c_s_square) * u_max;
                    } else {
                        int xmd = (nX + x - xDirection) % nX;
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


void stream_couette_arrays_opt1(int nX, int nY, int nZ, int direction_size, double c_s,
                                     double* previous_particle_distributions,
                                     double* particle_distributions,
                                     const int* directions,
                                     const double* weights,
                                     int* reverse_indexes) {
    double c_s_square = c_s * c_s;
    double inv_c_s_square = 1.0 / c_s_square;  // Pre-compute inverse of c_s_square
    double u_max = 0.1;  // Maximum velocity
    
    int nXY = nX * nY;
    int nXYZ = nXY * nZ;
    
    for (int x = 0; x < nX; x++) {
        for (int y = 0; y < nY; y++) {
            for (int z = 0; z < nZ; z++) {
                int baseIndex = x + y * nX + z * nXY;
                
                for (int i = 0; i < direction_size; i++) {
                    int index = baseIndex + i * nXYZ;
                    int reverseIndex = baseIndex + reverse_indexes[i] * nXYZ;
                    int xDirection = directions[3 * i];
                    int yDirection = directions[3 * i + 1];
                    int zDirection = directions[3 * i + 2];
                    
                    if (y == 0 && yDirection == 1) {
                        // Bottom Wall.
                        particle_distributions[index] = previous_particle_distributions[reverseIndex];
                    } else if (y == nY - 1 && yDirection == -1) {
                        // Top wall
                        particle_distributions[index] = previous_particle_distributions[reverseIndex] + xDirection * 2 * weights[i] * inv_c_s_square * u_max;
                    } else {
                        // Taylor-Green periodicity
                        int xmd = (nX + x - xDirection) % nX;
                        int ymd = y - yDirection;
                        int zmd = (nZ + z - zDirection) % nZ;
                        int otherIndex = xmd + ymd * nX + zmd * nXY + i * nXYZ;
                        particle_distributions[index] = previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}
