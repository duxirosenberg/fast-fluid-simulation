#include "LBM.h"
#include <immintrin.h> // For AVX instructions
#include <stddef.h>
#include <string.h>
#include <stdio.h>



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


void stream_couette_memcpy(struct LBMarrays* S) {
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
        int dir_nXYZ = i * nXYZ;                                //1 Intop
        int nX_min_xDir = nX - xDirection;                      //1 Intops
        int nZ_min_zDir = nZ - zDirection;                      //1 Intops

        int index = dir_nXYZ; 
        int reverseIndex = reverseIndex_nXYZ;
        int startY = 0;
        int endY = nY;

        if(yDirection == 1){//Bottom Wall, y=0
            for(int z = 0; z < nZ; z++) {           //Inner Loops; nZ*(2+nX* ((2/3 + 14/3)+(3/3+14/3))+((nY-2)*7))) Intops, nZ*nX*(1*(1/3)) flops
                int zmd = (nZ_min_zDir + z ) % nZ;                  //2 Intop
                memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[reverseIndex], nX * sizeof(double));
            }
            startY = 1;
        }else if (yDirection == -1){//TOP Wall, y=nY-1
            for(int z = 0; z < nZ; z++) {           //Inner Loops; nZ*(2+nX* ((2/3 + 14/3)+(3/3+14/3))+((nY-2)*7))) Intops, nZ*nX*(1*(1/3)) flops
                int zmd = (nZ_min_zDir + z ) % nZ;                  //2 Intop
                for(int x = 0; x < nX; x++) { 
                    int IIindex = x + (nY-1) * nX + z * nXY + i * nXYZ;
                    int RRreverseIndex = x + (nY-1) * nX + z * nXY + reverseIndex_nXYZ;
                    S->particle_distributions[IIindex] = S->previous_particle_distributions[RRreverseIndex] + temp; // 1 flop   // 1 FP READ, 1 FP WRITE
                }
            }
            endY = nY-1;
        }
        //Other cases
        if(xDirection == 0 && yDirection == 0 && zDirection == 0) {
            memcpy(&S->particle_distributions[dir_nXYZ], &S->previous_particle_distributions[dir_nXYZ], (sizeof(double)) * nXYZ);
        } else if(xDirection == 0 && yDirection == 0 && zDirection == -1) {
            memcpy(&S->particle_distributions[dir_nXYZ], &S->previous_particle_distributions[dir_nXYZ + nXY], (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dir_nXYZ + nXYZ - nXY], &S->previous_particle_distributions[dir_nXYZ], (sizeof(double)) * nXY);
        } else if(xDirection == 0 && yDirection == 0 && zDirection == 1) {
            memcpy(&S->particle_distributions[dir_nXYZ + S->nXY], &S->previous_particle_distributions[dir_nXYZ], (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dir_nXYZ], &S->previous_particle_distributions[dir_nXYZ + nXYZ - nXY], (sizeof(double)) * nXY);
        } else if(xDirection == 0 && yDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dir_nXYZ;
                int zIndex = z * nXY + dir_nXYZ;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + nX], (sizeof(double)) * (nXY - nX));
            }
        } else if(xDirection == 0 && yDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - yDirection) % nZ;
                int index = dir_nXYZ + nX;
                int otherZIndex = zmd * nXY + dir_nXYZ;
                memcpy(&S->particle_distributions[index + nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * (nXY - nX));
            }
        }else if(xDirection == -1) { 
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - yDirection) % nZ;
                int otherZIndex = zmd * nXY + dir_nXYZ;
                int zIndex = z * nXY + dir_nXYZ;
                for (int y = startY; y < endY; y++) { 
                    int ymd =  y - yDirection;
                    int otherIndex = ymd * nX + otherZIndex;
                    memcpy(&S->particle_distributions[dir_nXYZ], &S->previous_particle_distributions[otherIndex + 1], (sizeof(double)) * (nX - 1));
                    S->particle_distributions[dir_nXYZ + nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        }else if(xDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ_min_zDir + z) % nZ;
                int otherZIndex = zmd * nXY + dir_nXYZ;
                int zIndex = z * nXY + dir_nXYZ;
                for (int y = startY; y < endY; y++) {
                    int ymd = y - yDirection;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex + nX - 1];
                    memcpy(&S->particle_distributions[index + 1], &S->previous_particle_distributions[otherIndex], (sizeof(double)) * (nX - 1));
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
                }
            }
        }
    }
}








/*
void stream_couette_flattened_vec(struct LBMarrays* S){
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

    for(int i=0; i<q; i++){
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

        //startin vectorization
        int z = 0, y = 0, x = 0;
        __m256i z_vec = _mm256_setr_epi64x(0, 1, 2, 3);
        __m256i y_vec = _mm256_setr_epi64x(0, 1, 2, 3);
        __m256i x_vec = _mm256_setr_epi64x(0, 1, 2, 3);
        __m256i nX_min_xDir_vec = _mm256_set1_epi64x(nX - xDirection);
        __m256i nX_vec = _mm256_set1_epi64x(nX);
        
        int zmd = nZ_min_zDir % nZ;                   
        __m256i zmd_vec = _mm256_set1_epi64x(zmd);
        int ymd = -yDirection;                          //1 Intop
        __m256i ymd_vec = _mm256_set1_epi64x(ymd);
        int index0 = dir_nXYZ;                    
        __m256i indices = _mm256_setr_epi64x(dir_nXYZ, dir_nXYZ+1,dir_nXYZ+ 2,dir_nXYZ+ 3);       
        int reverseIndex0 = reverseIndex_nXYZ; 
        __m256i reverseIndices = _mm256_setr_epi64x(reverseIndex0, reverseIndex0+1,reverseIndex0+ 2,reverseIndex0+ 3);

        for(int j=0; j<nXYZ; j+=4){
            // xmd = (nX_min_xDir + x) % nX;
            //there is no AVX2 modulo: we can use the following trick: x % n = x - (x/n)*n
            __m256i x_sum_vec = _mm256_add_epi64(nX_min_xDir_vec, x_vec);
            //exact division
            
            __m256i x_mul_vec = _mm256_mullo_epi64(x_div_vec, nX_vec);
            __m256i xmd_vec = _mm256_sub_epi64(x_sum_vec, x_mul_vec);


        }
    }
}





void stream_couette_flattened_blaaaavx(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;
    double inv_c_s_square = 1.0 / c_s_square;
    double u_max = 0.1;
    double inv_c_s_square_u_max_times2 = inv_c_s_square * 0.1 * 2;
    int nX = S->nX;
    int nY = S->nY;
    int nZ = S->nZ;
    int q = S->direction_size;
    int nXY = S->nXY;
    int nXYZ = S->nXYZ;

    for (int i = 0; i < q; i++) {
        int xDirection = S->directions[3 * i];
        int yDirection = S->directions[3 * i + 1];
        int zDirection = S->directions[3 * i + 2];
        double temp = xDirection * S->weights[i] * inv_c_s_square_u_max_times2;
        int reverseIndex_nXYZ = S->reverse_indexes[i] * nXYZ;
        int yDirIS1 = yDirection == 1;
        int yDirISN1 = yDirection == -1;
        int dir_nXYZ = i * nXYZ;
        int nX_min_xDir = nX - xDirection;
        int nZ_min_zDir = nZ - zDirection;

        int z = 0, y = 0, x = 0;
        int zmd = nZ_min_zDir % nZ;
        int ymd = -yDirection;
        int index = dir_nXYZ;
        int reverseIndex = reverseIndex_nXYZ;

        // Load constants into SIMD registers
        __m256d temp_vec = _mm256_set1_pd(temp);
        __m256i nX_vec = _mm256_set1_epi32(nX);
        __m256i nXY_vec = _mm256_set1_epi32(nXY);
        __m256i nZ_vec = _mm256_set1_epi32(nZ);
        __m256i nX_min_xDir_vec = _mm256_set1_epi32(nX_min_xDir);
        __m256i yDirection_vec = _mm256_set1_epi32(yDirection);
        __m256i yDirIS1_vec = _mm256_set1_epi32(yDirIS1);
        __m256i yDirISN1_vec = _mm256_set1_epi32(yDirISN1);

        for (int j = 0; j < nXYZ; j += 4) { // Process 4 elements at a time
            //all input sizes are divisible by 4. 

            // Calculate indices
            __m256i x_vec = _mm256_setr_epi32(x, x + 1, x + 2, x + 3, 0, 0, 0, 0);
            __m256i x_sum_vec = _mm256_add_epi32(nX_min_xDir_vec, x_vec);
            __m256i x_div_vec = _mm256_srli_epi32(_mm256_mullo_epi32(x_sum_vec, _mm256_set1_epi32(0x55555555)), 31); // Approximate division
            __m256i x_mul_vec = _mm256_mullo_epi32(x_div_vec, nX_vec);
            __m256i xmd_vec = _mm256_sub_epi32(x_sum_vec, x_mul_vec);

            __m256i y_vec = _mm256_set1_epi32(y);
            __m256i zmd_vec = _mm256_set1_epi32(zmd);
            __m256i index_vec = _mm256_add_epi32(_mm256_set1_epi32(index), _mm256_setr_epi32(0, 1, 2, 3, 0, 0, 0, 0));
            __m256i reverseIndex_vec = _mm256_add_epi32(_mm256_set1_epi32(reverseIndex), _mm256_setr_epi32(0, 1, 2, 3, 0, 0, 0, 0));

            // Extract and gather particle distributions
            int reverseIndex_array[8];
            _mm256_storeu_si256((__m256i *)reverseIndex_array, reverseIndex_vec);
            __m256d particle_dist = _mm256_set_pd(
                S->previous_particle_distributions[reverseIndex_array[3]],
                S->previous_particle_distributions[reverseIndex_array[2]],
                S->previous_particle_distributions[reverseIndex_array[1]],
                S->previous_particle_distributions[reverseIndex_array[0]]
            );

            __m256d other_dist;

            // Apply conditions using masking
            __m256i cond1_mask = _mm256_cmpeq_epi32(y_vec, yDirIS1_vec);
            __m256i cond2_mask = _mm256_cmpeq_epi32(y_vec, _mm256_sub_epi32(_mm256_set1_epi32(nY), yDirISN1_vec));
            __m256d temp_masked = _mm256_and_pd(_mm256_castsi256_pd(cond2_mask), temp_vec);

            // Compute otherIndex
            __m256i otherIndex_vec = _mm256_add_epi32(
                _mm256_add_epi32(
                    _mm256_add_epi32(xmd_vec, _mm256_mullo_epi32(_mm256_set1_epi32(ymd), nX_vec)),
                    _mm256_mullo_epi32(zmd_vec, nXY_vec)),
                _mm256_set1_epi32(dir_nXYZ)
            );

            // Extract and gather other particle distributions
            int otherIndex_array[8];
            _mm256_storeu_si256((__m256i *)otherIndex_array, otherIndex_vec);
            other_dist = _mm256_set_pd(
                S->previous_particle_distributions[otherIndex_array[3]],
                S->previous_particle_distributions[otherIndex_array[2]],
                S->previous_particle_distributions[otherIndex_array[1]],
                S->previous_particle_distributions[otherIndex_array[0]]
            );

            // Blend results based on conditions
            __m256d result = _mm256_blendv_pd(other_dist, _mm256_add_pd(particle_dist, temp_masked), _mm256_castsi256_pd(cond1_mask));

            // Store results
            _mm256_storeu_pd(&S->particle_distributions[index + j], result);

            // Increment coordinates
            x += 4;
            index += 4;
            reverseIndex += 4;

            if (x >= nX) {
                x = 0;
                y++;
                ymd++;
                if (y >= nY) {
                    y = 0;
                    ymd = -yDirection;
                    z++;
                    zmd++;
                    if (zmd >= nZ) {
                        zmd = 0;
                    }
                }
            }
        }
    }
}


*/












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
