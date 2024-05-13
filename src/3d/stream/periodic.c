#include "LBM.h"

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

                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];
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