#include "LBM.h"
#include "string.h"
//#include <immintrin.h>

// Flops: 0
// Intops: xN * nY * nZ * q * 32
void stream_periodic_baseline(struct LBMarrays* S) {
    for(int x = 0; x < S->nX; x++) {
        for (int y = 0; y < S->nY; y++) {
            for (int z = 0; z < S->nZ; z++) {
                for (int i = 0; i < S->direction_size; i++) {
                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ; // 9 Intops
                    int xmd = (S->nX + x - S->directions[3 * i]) % S->nX; // 4 Intops
                    int ymd = (S->nY + y - S->directions[3 * i + 1]) % S->nY; // 5 Intops
                    int zmd = (S->nZ + z - S->directions[3 * i + 2]) % S->nZ; // 5 Intops
                    int otherIndex = xmd + ymd * S->nX + zmd * S->nX * S->nY + i * S->nX * S->nY * S->nZ; // 9 Intops

                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex] ;
                }
            }
        }
    }
}

// Flops: 0
// Intops: xN * nY * nZ * q * 32
void stream_periodic_arrays(int nX, int nY, int nZ, int direction_size,
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

///-------------OPTIMIZATION ONE -------------////


// Flops: 0
// Intops: nZ(1 + nY(1 + (nX * ( 2 + 19*q))))

void stream_periodic_O1(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int XY = SnX * SnY;
    int XYZ = XY * SnZ;
    int direction_size = S->direction_size;
    int index, idx, xmd, ymd, zmd;

    for (int z = 0; z < SnZ; z++) {
        int zSnXSnY = z * XY;
        for (int y = 0; y < SnY; y++) {
            int ySnX = y * SnX;
            for (int x = 0; x < SnX; x++) {
                int index = x + ySnX + zSnXSnY;
                for (int i = 0; i < direction_size; i++) {
                    idx = 3 * i;

                    xmd = (SnX + x - S->directions[idx]) % SnX; // 3 Intops
                    ymd = (SnY + y - S->directions[idx + 1]) % SnY; // 4 Intops
                    zmd = (SnZ + z - S->directions[idx + 2]) % SnZ; // 4 Intops

                    int otherIndex = xmd + ymd * SnX + zmd * XY + i * XYZ; //6 Intops

                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
                    index += XYZ;
                }
            }
        }
    }
}

///-------------OPTIMIZATION TWO--------------//


// Flops: 0
// Intops: 2 + q * (4 + (nZ * (4 + nY * (4 + 6 * nX))))



void stream_periodic_O2(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int XY = SnX * SnY;
    int XYZ = XY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* previous_particle_distributions = S->previous_particle_distributions;

    for (int i = 0; i < direction_size; i++) {
        int direction_x = directions[3 * i];
        int direction_y = directions[3 * i + 1];
        int direction_z = directions[3 * i + 2];

        int zmz = SnZ - direction_z;
        int ymy = SnY - direction_y;
        int xmx = SnX - direction_x;

        int iXYZ = i*XYZ;

        int index = iXYZ;
        

        for (int z = 0; z < SnZ; z++) {
            int zSnXSnY = z * XY;
            int zmd = (zmz + z) % SnZ;
            int zmdSnXSnY = zmd * XY;

            for (int y = 0; y < SnY; y++) {
                int ySnX = y * SnX;
                int ymd = (ymy + y) % SnY;
                int ymdSnX = ymd * SnX;

                for (int x = 0; x < SnX; x++) {
                    int xmd = (xmx + x) % SnX;
                    int otherIndex = xmd + ymdSnX + zmdSnXSnY + iXYZ;
                    particle_distributions[index] = previous_particle_distributions[otherIndex];
                    index++;
                }
            }
        }
    }
}

void stream_periodic_memcpy(struct LBMarrays* S, int time) {
    for (int i = 0; i < S->direction_size; i++) {
        // 6 * d
        int directionX = S->directions[3 * i];
        int directionY = S->directions[3 * i + 1];
        int directionZ = S->directions[3 * i + 2];
        int distIndex = i * S->nXYZ;
        // D3Q27: 13 iops
        // D3Q15: 13 iops
        // D2Q9: 1 iops
        if(directionX == 0 && directionY == 0 && directionZ == 0) {
            memcpy(&S->particle_distributions[distIndex], &S->previous_particle_distributions[distIndex], (sizeof(double)) * S->nXYZ);
        } else if(directionX == 0 && directionY == 0 && directionZ == -1) {
            memcpy(&S->particle_distributions[distIndex], &S->previous_particle_distributions[distIndex + S->nXY], (sizeof(double)) * (S->nXYZ - S->nXY));
            memcpy(&S->particle_distributions[distIndex + S->nXYZ - S->nXY], &S->previous_particle_distributions[distIndex], (sizeof(double)) * S->nXY);
        } else if(directionX == 0 && directionY == 0 && directionZ == 1) {
            memcpy(&S->particle_distributions[distIndex + S->nXY], &S->previous_particle_distributions[distIndex], (sizeof(double)) * (S->nXYZ - S->nXY));
            memcpy(&S->particle_distributions[distIndex], &S->previous_particle_distributions[distIndex + S->nXYZ - S->nXY], (sizeof(double)) * S->nXY);
        } else if(directionX == 0 && directionY == -1) {
            // D3Q27: z * 39 iops
            // D3Q15: z * 13 iops
            // D2Q9:  z * 13 iops
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + S->nX], (sizeof(double)) * (S->nXY - S->nX));
                memcpy(&S->particle_distributions[zIndex + S->nXY - S->nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * S->nX);
            }
        } else if(directionX == 0 && directionY == 1) {
            // D3Q27: z * 39 iops
            // D3Q15: z * 13 iops
            // D2Q9:  z * 13 iops
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                memcpy(&S->particle_distributions[zIndex + S->nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * (S->nXY - S->nX));
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + S->nXY - S->nX], (sizeof(double)) * S->nX);
            }
        } else if(directionX == -1) {
            // D3Q27: 9 * z * (7 + 12 * y) iops
            // D3Q15: 5 * z * (7 + 12 * y) iops
            // D2Q9:  3 * z * (7 + 12 * y) iops
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                for (int y = 0; y < S->nY; y++) {
                    int ymd = (S->nY + y - directionY) % S->nY;
                    int otherIndex = ymd * S->nX + otherZIndex;
                    int index = y * S->nX + zIndex;
                    memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[otherIndex + 1], (sizeof(double)) * (S->nX - 1));
                    S->particle_distributions[index + S->nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        } else if(directionX == 1) {
            // D3Q27: 9 * z * (7 + 12 * y) iops
            // D3Q15: 5 * z * (7 + 12 * y) iops
            // D2Q9:  3 * z * (7 + 12 * y) iops
            for (int z = 0; z < S->nZ; z++) {
                int zmd = (S->nZ + z - directionZ) % S->nZ;
                int otherZIndex = zmd * S->nXY + distIndex;
                int zIndex = z * S->nXY + distIndex;
                for (int y = 0; y < S->nY; y++) {
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



/*void stream_periodic_loopunrolling(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* previous_particle_distributions = S->previous_particle_distributions;


    for (int z = 0; z < SnZ; z++) {
        for (int y = 0; y < SnY; y++) {
            for (int x = 0; x < SnX; x++) {
                int base_index = x + y * SnX + z * SnXSnY;

                int i;
                for (i = 0; i <= direction_size - 4; i += 4) {
                    int idx0 = 3 * i;
                    int idx1 = 3 * (i + 1);
                    int idx2 = 3 * (i + 2);
                    int idx3 = 3 * (i + 3);

                    int xmd0 = (SnX + x - directions[idx0]) % SnX;
                    int ymd0 = (SnY + y - directions[idx0 + 1]) % SnY;
                    int zmd0 = (SnZ + z - directions[idx0 + 2]) % SnZ;
                    int xmd1 = (SnX + x - directions[idx1]) % SnX;
                    int ymd1 = (SnY + y - directions[idx1 + 1]) % SnY;
                    int zmd1 = (SnZ + z - directions[idx1 + 2]) % SnZ;
                    int xmd2 = (SnX + x - directions[idx2]) % SnX;
                    int ymd2 = (SnY + y - directions[idx2 + 1]) % SnY;
                    int zmd2 = (SnZ + z - directions[idx2 + 2]) % SnZ;
                    int xmd3 = (SnX + x - directions[idx3]) % SnX;
                    int ymd3 = (SnY + y - directions[idx3 + 1]) % SnY;
                    int zmd3 = (SnZ + z - directions[idx3 + 2]) % SnZ;

                    int otherIndex0 = xmd0 + ymd0 * SnX + zmd0 * SnXSnY + i * SnXSnYSnZ;
                    int otherIndex1 = xmd1 + ymd1 * SnX + zmd1 * SnXSnY + (i + 1) * SnXSnYSnZ;
                    int otherIndex2 = xmd2 + ymd2 * SnX + zmd2 * SnXSnY + (i + 2) * SnXSnYSnZ;
                    int otherIndex3 = xmd3 + ymd3 * SnX + zmd3 * SnXSnY + (i + 3) * SnXSnYSnZ;

                    particle_distributions[base_index + i * SnXSnYSnZ] = previous_particle_distributions[otherIndex0];
                    particle_distributions[base_index + (i + 1) * SnXSnYSnZ] = previous_particle_distributions[otherIndex1];
                    particle_distributions[base_index + (i + 2) * SnXSnYSnZ] = previous_particle_distributions[otherIndex2];
                    particle_distributions[base_index + (i + 3) * SnXSnYSnZ] = previous_particle_distributions[otherIndex3];
                }

                // Handle any remaining elements
                for (; i < direction_size; i++) {
                    int idx = 3 * i;
                    int xmd = (SnX + x - directions[idx]) % SnX;
                    int ymd = (SnY + y - directions[idx + 1]) % SnY;
                    int zmd = (SnZ + z - directions[idx + 2]) % SnZ;

                    int otherIndex = xmd + ymd * SnX + zmd * SnXSnY + i * SnXSnYSnZ;
                    particle_distributions[base_index + i * SnXSnYSnZ] = previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}*/







/*void stream_periodic_O3(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int SnXSnY = SnX * SnY;
    int SnXSnYSnZ = SnXSnY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* previous_particle_distributions = S->previous_particle_distributions;
    
    for (int z = 0; z < SnZ; z++) {
        for (int y = 0; y < SnY; y++) {
            for (int x = 0; x < SnX; x++) {
                int base_index = x + y * SnX + z * SnXSnY;

                for (int i = 0; i < direction_size; i += 4) {
                    __m256i idx = _mm256_loadu_si256((__m256i*)&directions[3 * i]);
                    __m256i x_off = _mm256_and_si256(_mm256_add_epi32(_mm256_set1_epi32(x), _mm256_sub_epi32(_mm256_set1_epi32(SnX), _mm256_set1_epi32(directions[3 * i]))), _mm256_set1_epi32(SnX - 1));
                    __m256i y_off = _mm256_and_si256(_mm256_add_epi32(_mm256_set1_epi32(y), _mm256_sub_epi32(_mm256_set1_epi32(SnY), _mm256_set1_epi32(directions[3 * i + 1]))), _mm256_set1_epi32(SnY - 1));
                    __m256i z_off = _mm256_and_si256(_mm256_add_epi32(_mm256_set1_epi32(z), _mm256_sub_epi32(_mm256_set1_epi32(SnZ), _mm256_set1_epi32(directions[3 * i + 2]))), _mm256_set1_epi32(SnZ - 1));

                    __m256i other_index = _mm256_add_epi32(_mm256_add_epi32(_mm256_mullo_epi32(z_off, _mm256_set1_epi32(SnXSnY)), _mm256_mullo_epi32(y_off, _mm256_set1_epi32(SnX))), x_off);

                    for (int j = 0; j < 4; ++j) {
                        int idx = _mm256_extract_epi32(other_index, j) + (i + j) * SnXSnYSnZ;
                        particle_distributions[base_index + (i + j) * SnXSnYSnZ] = previous_particle_distributions[idx];
                    }
                }
            }
        }
    }
}
*/




/*void stream_periodic_ZYXI(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int XY = SnX * SnY;
    int XYZ = XY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* previous_particle_distributions = S->previous_particle_distributions;

 
    for (int z = 0; z < SnZ; z++) {
        int zSnXSnY = z * XY;
        for (int y = 0; y < SnY; y++) {
            int ySnX = y * SnX;
            for (int x = 0; x < SnX; x++) {
                int index = x + y * SnX + z * XY;
                for (int i = 0; i < direction_size; i++) {
                    int idx = 3 * i;
                    int xmd = (SnX + x - directions[idx]) % SnX;
                    int ymd = (SnY + y - directions[idx + 1]) % SnY;
                    int zmd = (SnZ + z - directions[idx + 2]) % SnZ;

                    int otherIndex = xmd + ymd * SnX + zmd * XY + i * XYZ;
                    particle_distributions[index ] = previous_particle_distributions[otherIndex];
                    index += XYZ;
                }
            }
        }
    }
}*/

/*void stream_periodic_O21(struct LBMarrays* S) {
    int SnX = S->nX;
    int SnY = S->nY;
    int SnZ = S->nZ;
    int XY = SnX * SnY;
    int XYZ = XY * SnZ;
    int direction_size = S->direction_size;
    const int* directions = S->directions;
    double* particle_distributions = S->particle_distributions;
    double* previous_particle_distributions = S->previous_particle_distributions;

    for (int i = 0; i < direction_size; i++) {
        int direction_x = directions[3 * i];
        int direction_y = directions[3 * i + 1];
        int direction_z = directions[3 * i + 2];

        int zmz = SnZ - direction_z;
        int ymy = SnY - direction_y;
        int xmx = SnX - direction_x;

        for (int z = 0; z < SnZ; z++) {
            int zSnXSnY = z * XY;
            int zmd = (zmz + z) % SnZ;
            int zmdSnXSnY = zmd * XY;

            for (int y = 0; y < SnY; y++) {
                int ySnX = y * SnX;
                int ymd = (ymy + y) % SnY;
                int ymdSnX = ymd * SnX;

                for (int x = 0; x < SnX; x++) {
                    int index = x + ySnX + zSnXSnY + i * XYZ;
                    int xmd = (xmx + x) % SnX;
                    int otherIndex = xmd + ymdSnX + zmdSnXSnY + i * XYZ;
                    particle_distributions[index] = previous_particle_distributions[otherIndex];
                }
            }
        }
    }
}*/