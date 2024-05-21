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
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == S->nY - 1 && yDirection == -1) { // 5 Flops (nX * nZ * direction_size / 3) times
                        double u_max = 0.1;
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection * 2 * S->weights[i] / (c_s_square) * u_max;
                    } else { // 16 Intops, (nX*nZ*( (nY-2) + 4*direction_size/3 )
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
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == nY - 1 && yDirection == -1) { // 4 Intops, 3 Flops (nX * nZ * direction_size / 3) times
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection *  S->weights[i] * inv_c_s_square_u_max_times2;
                    } else { // 16 Intops, 1 FP READ, 1 FP WRITE (nX * nZ * ((nY-2) + (4*direction_size / 3))) times
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




// Flops: 4 + q(2 + nZ*nX*(1*(1/3)))
// Intops: q*(11+nZ*(2+nX* ((2/3 + 14/3)+(3/3+14/3))+((nY-2)*7))))
// INTS READ: q*4
// INTS Written: 0
// FLOPS READ: q*(nZ*nY*NX + 1)
// FLOPS WRITTEN: q(nZ*nY*NX)
// Loop Unrolling: no effect, no reduction, no dependencies, compiler can do it
// getting if/else outside of loop: seems to do worse, compiler does something smarter.  
void stream_couette_loop_structure(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;                        //1 Flop
    double inv_c_s_square = 1.0 / c_s_square;                   //1 Flop
    double u_max = 0.1;  // Maximum velocity
    double inv_c_s_square_u_max_times2 = inv_c_s_square * 0.1*2; //2 Flops
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    int nXY = S->nXY;
    int nXYZ = S->nXYZ;

    for(int i=0;i<q; i++){   //OUTER LOOP: q iterations
        int xDirection = S->directions[3 * i];                  //1 Intop,      1 INT READ
        int yDirection = S->directions[3 * i + 1];              //2 Intops,     1 INT READ
        int zDirection = S->directions[3 * i + 2];              //2 Intops,     1 INT READ
        double temp = xDirection * S->weights[i] * inv_c_s_square_u_max_times2; //2 Flops   1 FL READ
        int reverseIndex_nXYZ = S->reverse_indexes[i] * nXYZ;   //1 Intops      1 INT READ
        int yDirIS1 = yDirection == 1;                          //1 Intop
        int yDirISN1 = yDirection == -1;                        //1 Intop
        int dir_nXYZ = i * nXYZ;                                //1 Intop
        int nX_min_xDir = nX - xDirection;                      //1 Intops
        int nZ_min_zDir = nZ - zDirection;                      //1 Intops

        int index = dir_nXYZ; 
        int reverseIndex = reverseIndex_nXYZ;

        for(int z = 0; z < nZ; z++) {           //Inner Loops; nZ*(2+nX* ((2/3 + 14/3)+(3/3+14/3))+((nY-2)*7))) Intops, nZ*nX*(1*(1/3)) flops
            int zmd = (nZ_min_zDir + z ) % nZ;                  //2 Intop
            for (int y = 0; y < nY; y++) {                      
                int ymd = y - yDirection;                       //1 Intop   
                for (int x = 0; x < nX; x++) {
                    if (y == 0 && yDirIS1) {                    //2 Intops   (yDIRIS1 is true in 1/3 of directions)
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];          // 1 FP READ, 1 FP WRITE
                    } else if (y == nY - 1 && yDirISN1) {       // 3 Intops  (yDirISN1 is 1/3 in 1/3 of directions)
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + temp; // 1 flop   // 1 FP READ, 1 FP WRITE
                    } else {                                   
                        int xmd = (nX_min_xDir + x ) % nX;      //2 Intops
                        int otherIndex = xmd + ymd * nX + zmd * nXY + dir_nXYZ;     //5 Intops
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];      // 1 FP READ, 1 FP WRITE
                    }
                    index++;                                     //1 Intop
                    reverseIndex++;                              //1 Intop
                }
            }
        }
    }
}









//Flops: q * (nX * nY * nZ + 2) + 4
//Intops: q * (12 +(nZ*(2+nX* ((2/3 + 10/3)+(3/3+10/3))+((nY-2)*5)))) +  (nX*nY*nZ + 4*nZ*nY + 5* NZ + 1))
// INTS READ: q*4
// INTS Written: 0
// FLOPS READ: q*(nZ*nY*NX + 1)
// FLOPS WRITTEN: q(nZ*nY*NX)
// Loop Unrolling: no effect, no reduction, no dependencies, compiler can do it
// getting if/else outside of loop: seems to do worse, compiler does something smarter.  
void stream_couette_flattened(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;                //1 Flop
    double inv_c_s_square = 1.0 / c_s_square;           //1 Flop   
    double u_max = 0.1;  
    double inv_c_s_square_u_max_times2 = inv_c_s_square * 0.1 * 2;  //2 Flops
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    int nXY = S->nXY;
    int nXYZ = S->nXYZ;

    for (int i = 0; i < q; i++) {                       //OUTER LOOP: q iterations (13 + 26*nX*nY*nZ* Intops,  2 + nX*nY*nZ Flops)               
        int xDirection = S->directions[3 * i];          //1 Intop       // 1 INT READ
        int yDirection = S->directions[3 * i + 1];      //2 Intops      // 1 INT READ
        int zDirection = S->directions[3 * i + 2];      //2 Intops      // 1 INT READ
        double temp = xDirection * S->weights[i] * inv_c_s_square_u_max_times2; //2 Flops   // 1 FP READ
        int reverseIndex_nXYZ = S->reverse_indexes[i] * nXYZ;   //1 Intops                  // 1 INT READ
        int yDirIS1 = yDirection == 1;                  //1 Intop
        int yDirISN1 = yDirection == -1;                //1 Intop
        int dir_nXYZ = i * nXYZ;                        //1 Intop
        int nX_min_xDir = nX - xDirection;              //1 Intops
        int nZ_min_zDir = nZ - zDirection;              //1 Intops

        int z = 0, y = 0, x = 0;
        int zmd = nZ_min_zDir % nZ;                     //1 Intop
        int ymd = -yDirection;                          //1 Intop
        int index = dir_nXYZ;                           
        int reverseIndex = reverseIndex_nXYZ;          

        for (int j = 0; j < nXYZ; j++) {                //Inner Loops; nZ*(2+nX* ((2/3 + 10/3)+(3/3+10/3))+((nY-2)*5))) Intops, nZ*nX*(1*(1/3)) flops
            int xmd = (nX_min_xDir + x) % nX;           //2 Intops

            if (y == 0 && yDirIS1) {                    //2 Intops 
                S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];      // 1 FP READ, 1 FP WRITE
            } else if (y == nY - 1 && yDirISN1) {       //3 Intops
                S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + temp; //1 Flop // 1 FP READ, 1 FP WRITE
            } else {
                int otherIndex = xmd + ymd * nX + zmd * nXY + dir_nXYZ; //5 Intops
                S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];  // 1 FP READ, 1 FP WRITE
            }

            // Increment coordinates
            x++;                                        //1 Intop
            index++;                                    //1 Intop
            reverseIndex++;                             //1 Intop

            //this below has a total of nX*nY*nZ + 4*nZ*nY + 5* NZ + 1 Intops
            if (x == nX) {                              //1 Intop
                //true only nZ*nY times
                x = 0;                                  //1 Intop
                y++;                                    //1 Intop
                ymd++;                                  //1 Intop
                if (y == nY) {                          //1 Intop
                    //true only nZ times
                    y = 0;                              //1 Intop
                    ymd = -yDirection;                  //1 Intop
                    z++;                                //1 Intop
                    zmd++;                              //1 Intop
                    if (zmd == nZ) {                    //1 Intop
                        //true only once
                        zmd = 0;                        //1 Intop
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
