#include "LBM.h"
#include <immintrin.h> // For AVX instructions

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







// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
void stream_couette_code_motion(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;
    double inv_c_s_square = 1.0 / c_s_square;  
    double u_max = 0.1;  // Maximum velocity
    double inv_c_s_square_u_max_times2 = inv_c_s_square * 0.1*2;
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    int nXY = S->nXY;
    int nXYZ = S->nXYZ;

    for(int x = 0; x < nX; x++) {
        for (int y = 0; y < nY; y++) {
            for (int z = 0; z < nZ; z++) {
                for (int i = 0; i < q; i++) {
                    int index = x + y * nX + z * nXY + i * nXYZ; // 9 Intops
                    int reverseIndex =  x + y * nX + z * nXY + S->reverse_indexes[i] * nXYZ; // 9 Intops
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
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection *  S->weights[i] * inv_c_s_square_u_max_times2;
                    } else { // 16 Intops, 0 Flops
                        // Chapter 13 Taylor-Green periodicity from same book.
                        int xmd = (nX + x - xDirection) % nX;
                        int ymd = y - yDirection;
                        int zmd = (nZ + z - zDirection) % nZ;
                        int otherIndex = xmd + ymd * nX + zmd * nX * nY + i * nX * nY * nZ;
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    }
                }
            }
        }
    }
}




// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
void stream_couette_loop_structure(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;
    double inv_c_s_square = 1.0 / c_s_square;  
    double u_max = 0.1;  // Maximum velocity
    double inv_c_s_square_u_max_times2 = inv_c_s_square * 0.1*2;
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    int nXY = S->nXY;
    int nXYZ = S->nXYZ;

    for(int i=0;i<q; i++){
        int xDirection = S->directions[3 * i];
        int yDirection = S->directions[3 * i + 1];
        int zDirection = S->directions[3 * i + 2];
        double temp = xDirection * S->weights[i] * inv_c_s_square_u_max_times2;
        int reverseIndex_nXYZ = S->reverse_indexes[i] * nXYZ;
        int yDirIS1 = yDirection == 1;
        int yDirISN1 = yDirection == -1;
        int dir_nXYZ = i * nXYZ;

        int index = dir_nXYZ; 
        int reverseIndex = reverseIndex_nXYZ;

        int nX_min_xDir = nX - xDirection;
        int nZ_min_zDir = nZ - zDirection;
        
        int xmd = nX_min_xDir;

        for(int z = 0; z < nZ; z++) {
            int zmd = (nZ_min_zDir + z ) % nZ;
            for (int y = 0; y < nY; y++) {
                int ymd = y - yDirection;
                for (int x = 0; x < nX; x++) {
                    if (y == 0 && yDirIS1) { //2 Intops, 1 Memory Op
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == nY - 1 && yDirISN1) { // 3 Intops, 1 Memory Op
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + temp;
                    } else { // 16 Intops, 0 Flops
                        int xmd = (nX_min_xDir + x ) % nX;
                        int otherIndex = xmd + ymd * nX + zmd * nXY + dir_nXYZ;
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    }
                    index++;
                    reverseIndex++;
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
