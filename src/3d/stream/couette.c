#include "LBM.h"
#include <immintrin.h> // For AVX instructions
#include <stddef.h>
#include <string.h>



/*BASELINE:
Flops: 5 * (nX * nZ * q / 3)
Intops: nX * nZ * q * (39 * nY + 32 / 3)
*/
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


/* OPTIMIZATION 1: CODE MOTION

- I moved the calculation of as many constants out of the loop as possible without changing the loop order.
- This optimization had a minimal effect on runtime, I suppose compiler does this well. 

OP COUNTING: 
IntOps: q* [ 12nXYZ + nXYZ-2(nX*nZ*(1/3))]
Flops: 4 + q*(3*nXZ*(1/3))
Doubles Read: 2*q*nXYZ + 27(for the weights)
Ints Read: 6 + 3*q(directions)
Doubles Written: q*nXYZ
*/
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
                    int index = x + y * nX + z * nXY + i * nXYZ;                                // 6 Intops
                    int reverseIndex =  x + y * nX + z * nXY + S->reverse_indexes[i] * nXYZ;    // 6 Intops
                    int xDirection = S->directions[3 * i];
                    int yDirection = S->directions[3 * i + 1];
                    int zDirection = S->directions[3 * i + 2];
                    if (y == 0 && yDirection == 1) {// Bottom Wal: (nXZ * q / 3) times
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];
                    } else if (y == nY - 1 && yDirection == -1) { // 3 Flops (nXZ * q / 3) times
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + xDirection *  S->weights[i] * inv_c_s_square_u_max_times2;
                    } else { // 13 Intops nXYZ-2(nX*nZ*(1/3)q) times
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


/*OPTIMIZATION 2: LOOP STRUCTURE

- rearranged the loops in a smarter way, while keeping the logic the same
- Loop unrolling had no effect. I think this is because there are no reductions or dependencies in the loop, 
hence the compiler can do it perfeclty well, maybe even better.
- taking the if/else out of the loop also had no effect, the compiler seems to do something smarter.

OP COUNTING: 
IntOps: q * (11 + [nZ*(2+ nY*(1 +nX( (1/3)*7 + 2)))] )
Flops: 4 + q * (2 + nZ*nX*(1/3))
Doubles Read: q * (1 + 2*nZ*nY*nX) 
Ints Read: q * 4
Doubles Written: q * nZ*nY*nX
*/  
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

        /* OP COUNTING IN THE Below Loops: 

        Intops: nZ*(2+ nY*(1 +nX( (1/3)*7 + 2)))
        Flops: nZ*nY*NX*(1/3)
        Doubles Read: 2*nZ*nY*nX
        Doubles Written: nZ*nY*nX
        */
        for(int z = 0; z < nZ; z++) {           
            int zmd = (nZ_min_zDir + z ) % nZ;                  //2 Intop
            for (int y = 0; y < nY; y++) {                      
                int ymd = y - yDirection;                       //1 Intop   
                for (int x = 0; x < nX; x++) {
                    if (y == 0 && yDirIS1) {                //(1/3)* 2 double reads, 1 double write                
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex];         
                    } else if (y == nY - 1 && yDirISN1) {   //(1/3)* 2 double reads, 1 double write, 1 flop
                        S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + temp; 
                    } else {                                //(1/3) * 7 intops, 2 double read, 1 double write               
                        int xmd = (nX_min_xDir + x ) % nX;      
                        int otherIndex = xmd + ymd * nX + zmd * nXY + dir_nXYZ;    
                        S->particle_distributions[index] = S->previous_particle_distributions[otherIndex];     
                    }
                    index++;                                     //1 Intop
                    reverseIndex++;                              //1 Intop
                }
            }
        }
    }
}


/*OPTIMIZATION 3: USING MEMCPY WHEREVER POSSIBLE

IntOps: 9+ q * [(1/3) * (nZ*6 Intops) + (1/3)*(10+nX/2 Intops) + 4 + 4]
Flops:  6 + q * [(1/3)*nZ*nX*Flops]
Doubles Read: 1 + q *[2*(1/3)*nZ*nX + 2*(1/3)*nZ*nX + 6*(1/q)*nXYZ + 4 *(1/15)_or_(3/27)*nZ*(nXY-nX) + 2*(1/3)*nZ*(nY-(1/3))*nX]
Ints Read: 10 
Doubles Written: q * [(1/3)*nZ*nX + (1/3)*nZ*nX + 3*(1/q)*nXYZ + 2*(1/15)_or_(3/27)*nZ*(nXY-nX) + *(1/3)*nZ*(nY-(1/3))*nX]
*/
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
        int reverseIndex_nXYZ = S->reverse_indexes[i] * nXYZ;   //1 Intops      1 INT READ
        int dirIndex = i * nXYZ;                                //1 Intop
        int nX_min_xDir = nX - xDirection;                      //1 Intops
        int nZ_min_zDir = nZ - zDirection;                      //1 Intops

        int startY = 0;
        int endY = nY;

        if(yDirection == 1){//Bottom Wall, 
            //y=0
            for(int z = 0; z < nZ; z++) {           //Inner Loops; nZ*(2+nX* ((2/3 + 14/3)+(3/3+14/3))+((nY-2)*7))) Intops, nZ*nX*(1*(1/3)) flops
                int zmd = (nZ_min_zDir + z ) % nZ;                  //2 Intop
                int index = z * nXY + dirIndex; // 9 Intops
                int reverseIndex =  z * nXY + reverseIndex_nXYZ; // 9 Intops
                memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[reverseIndex], nX * sizeof(double));
            }
            startY++;
        }else if (yDirection == -1){//TOP Wall,
            double temp = xDirection * S->weights[i] * inv_c_s_square_u_max_times2; //2 Flops   1 FL READ
            //y=nY-1
            for(int z = 0; z < nZ; z++) {           //Inner Loops; nZ*(2+nX* ((2/3 + 14/3)+(3/3+14/3))+((nY-2)*7))) Intops, nZ*nX*(1*(1/3)) flops
                int zIndex = (nY-1) * nX + z * nXY + dirIndex;
                int zReverseIndex = (nY-1) * nX + z * nXY + reverseIndex_nXYZ;
                for(int x = 0; x < nX; x++) { 
                    int index = zIndex + x;
                    int reverseIndex = zReverseIndex + x;
                    S->particle_distributions[index] = S->previous_particle_distributions[reverseIndex] + temp; // 1 flop   // 1 FP READ, 1 FP WRITE
                }
            }
            endY--;
        }
        //other cases
        if (xDirection == 0 && yDirection == 0 && zDirection == 0) {
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex], (sizeof(double)) * nXYZ);
        } else if (xDirection == 0 && yDirection == 0 && zDirection == -1) {
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex + nXY], (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dirIndex + nXYZ - nXY], &S->previous_particle_distributions[dirIndex], (sizeof(double)) * nXY);
        } else if (xDirection == 0 && yDirection == 0 && zDirection == 1) {
            memcpy(&S->particle_distributions[dirIndex + nXY], &S->previous_particle_distributions[dirIndex], (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex + nXYZ - nXY], (sizeof(double)) * nXY);
        } else if (xDirection == 0 && yDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ_min_zDir + z) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + nX], (sizeof(double)) * (nXY - nX));
            }
        } else if (xDirection == 0 && yDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                memcpy(&S->particle_distributions[zIndex + nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * (nXY - nX));
            }
        } else if (xDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                for (int y = startY; y < endY; y++) { 
                    int ymd =  y - yDirection;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[otherIndex + 1], (sizeof(double)) * (nX - 1));
                    S->particle_distributions[index + nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        } else if (xDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = y - yDirection;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    memcpy(&S->particle_distributions[index + 1], &S->previous_particle_distributions[otherIndex], (sizeof(double)) * (nX - 1));
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex + nX - 1];
                }
            }
        }
    }
}


/*OPTIMIZATION 4: USING AVX INSTRUCITONS

IntOps: 9+ q * [(1/3) * (nZ*6 Intops) + (1/3)*(10+nX/2 Intops) + 4 + 4]
Flops:  6 + q * [(1/3)*nZ*nX*Flops]
Doubles Read: 1 + q *[2*(1/3)*nZ*nX + 2*(1/3)*nZ*nX + 6*(1/q)*nXYZ + 4 *(1/15)_or_(3/27)*nZ*(nXY-nX) + 2*(1/3)*nZ*(nY-(1/3))*nX]
Ints Read: 10 
Doubles Written: q * [(1/3)*nZ*nX + (1/3)*nZ*nX + 3*(1/q)*nXYZ + 2*(1/15)_or_(3/27)*nZ*(nXY-nX) + *(1/3)*nZ*(nY-(1/3))*nX]
*/
void stream_couette_avx(struct LBMarrays* S) {
    double c_s_square = S->c_s * S->c_s;                        //1 Flop
    double inv_c_s_square = 1.0 / c_s_square;                   //1 Flop
    double u_max = 0.1;  // Maximum velocity
    double inv_c_s_square_u_max_times2 = inv_c_s_square * 0.1*2; //2 Flops
    int nX = S->nX;                                              // Read 1 int
    int nY = S->nY;                                              // Read 1 int
    int nZ = S->nZ;                                              // Read 1 int
    int q = S->direction_size;                                   // Read 1 int
    int nXY = S->nXY;                                            // Read 1 int
    int nXYZ = S->nXYZ;                                          // Read 1 int

    for(int i=0;i<q; i++){   //OUTER LOOP: q iterations
        int xDirection = S->directions[3 * i];                  //1 Intop,      Read 1 int
        int yDirection = S->directions[3 * i + 1];              //2 Intops,     Read 1 int
        int zDirection = S->directions[3 * i + 2];              //2 Intops,     Read 1 int
        double temp = xDirection * S->weights[i] * inv_c_s_square_u_max_times2; //2 Flops  Read 1 double
        int reverseIndex_nXYZ = S->reverse_indexes[i] * nXYZ;   //1 Intops      Read 1 int
        int dirIndex = i * nXYZ;                                //1 Intop
        int nX_min_xDir = nX - xDirection;                      //1 Intops
        int nZ_min_zDir = nZ - zDirection;                      //1 Intops

        int startY = 0;
        int endY = nY;

        /*
        OPS COUNTING IN THE BELOW IF/ELSE BLOCKS:
        Case yDir=1:    (1/3) * (nZ*6 Intops), 2*(1/3)*nZ*nX*sizeof(double) Bytes Read, (1/3)*nZ*nX*sizeof(double) Bytes Written
        Case yDir=-1:   (1/3)*(10+nX/2 Intops),  (1/3)*nZ*nX*Flops,  2*(1/3)*nZ*nX*sizeof(double) Bytes Read,  (1/3)*nZ*nX*sizeof(double) Bytes Written
        Case (0,0,0):   2*(1/q)*nXYZ*sizeof(double) Bytes Read, (1/q)*nXYZ*sizeof(double) Bytes Written
        Case (0,0,-1):  2*(1/q)*nXYZ*sizeof(double) Bytes Read, (1/q)*nXYZ*sizeof(double) Bytes Written, (1/q)*4 Intops
        Case (0,0,1):   2*(1/q)*nXYZ*sizeof(double) Bytes Read, (1/q)*nXYZ*sizeof(double) Bytes Written, (1/q)*4 Intops
        Case (0,-1,*):   2*(1/15)_or_(3/27)*nZ*(nXY-nX)*sizeof(double) Bytes Read, (1/15)_or_(3/27)*nZ*(nXY-nX)* sizeof(double) Bytes Written, (1/15)_or_(3/27)*nZ*5 Intops
        Case (0,1,*):    2*(1/15)_or_(3/27)*nZ*(nXY-nX)*sizeof(double) Bytes Read, (1/15)_or_(3/27)*nZ*(nXY-nX)* sizeof(double) Bytes Written, (1/15)_or_(3/27)*nZ*5 Intops
        Case (-1,*,*):     (1/3)*nZ*(5+nY*5) Intops, 2*(1/3)*nZ*(nY-(1/3))*nX*sizeof(double) Bytes Read, (1/3)*nZ*(nY-(1/3))*nX*sizeof(double) Bytes Written
        Case (1,*,*):      (1/3)*nZ*(5+nY*5) Intops, 2* (1/3)*nZ*(nY-(1/3))*nX*sizeof(double) Bytes Read, (1/3)*nZ*(nY-(1/3))*nX*sizeof(double) Bytes Written
        
        TOTAL:
        IntOps: (1/3) * (nZ*6 Intops) + (1/3)*(10+nX/2 Intops) + 4 + 4
        Flops:  (1/3)*nZ*nX*Flops
        Doubles Read: 2*(1/3)*nZ*nX + 2*(1/3)*nZ*nX + 6*(1/q)*nXYZ + 4 *(1/15)_or_(3/27)*nZ*(nXY-nX) + 2*(1/3)*nZ*(nY-(1/3))*nX
        Bytes Written: (1/3)*nZ*nX + (1/3)*nZ*nX + 3*(1/q)*nXYZ + 2*(1/15)_or_(3/27)*nZ*(nXY-nX) + *(1/3)*nZ*(nY-(1/3))*nX
        */


        if(yDirection == 1){//BOTTOM WALL,  
            //y=0
            for(int z = 0; z < nZ; z++) {          
                int zmd = (nZ_min_zDir + z ) % nZ;                  
                int index = z * nXY + dirIndex;                     
                int reverseIndex =  z * nXY + reverseIndex_nXYZ;   
                memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[reverseIndex], nX * sizeof(double));
            }
            startY++;
        }else if (yDirection == -1){//TOP Wall, 
            //y=nY-1
            for(int z = 0; z < nZ; z++) {      
                int zIndex = (nY-1) * nX + z * nXY + dirIndex;                      
                int zReverseIndex = (nY-1) * nX + z * nXY + reverseIndex_nXYZ;      
            
                __m256d temp_vec = _mm256_set1_pd(temp); // Load temp into an AVX register
                for(int x = 0; x < nX; x+=4) {                                      
                    int index = zIndex + x;                                         
                    int reverseIndex = zReverseIndex + x;                           
                    // Load 4 floats from previous_particle_distributions (reverseIndex until reverseIndex + 3) 
                    __m256d prev_part_dist_vec = _mm256_loadu_pd(&S->previous_particle_distributions[reverseIndex]);
                    __m256d result_vec = _mm256_add_pd(prev_part_dist_vec, temp_vec);
                    // Store the result back into particle_distributions
                    _mm256_storeu_pd(&S->particle_distributions[index], result_vec);                                   
                }
            }
            endY--;
        }
        //other cases
        if (xDirection == 0 && yDirection == 0 && zDirection == 0) {
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex], (sizeof(double)) * nXYZ);
        } else if (xDirection == 0 && yDirection == 0 && zDirection == -1) {
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex + nXY], (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dirIndex + nXYZ - nXY], &S->previous_particle_distributions[dirIndex], (sizeof(double)) * nXY);
        } else if (xDirection == 0 && yDirection == 0 && zDirection == 1) {
            memcpy(&S->particle_distributions[dirIndex + nXY], &S->previous_particle_distributions[dirIndex], (sizeof(double)) * (nXYZ - nXY));
            memcpy(&S->particle_distributions[dirIndex], &S->previous_particle_distributions[dirIndex + nXYZ - nXY], (sizeof(double)) * nXY);
        } else if (xDirection == 0 && yDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ_min_zDir + z) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                memcpy(&S->particle_distributions[zIndex], &S->previous_particle_distributions[otherZIndex + nX], (sizeof(double)) * (nXY - nX));
            }
        } else if (xDirection == 0 && yDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                memcpy(&S->particle_distributions[zIndex + nX], &S->previous_particle_distributions[otherZIndex], (sizeof(double)) * (nXY - nX));
            }
        } else if (xDirection == -1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd =  y - yDirection;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    memcpy(&S->particle_distributions[index], &S->previous_particle_distributions[otherIndex + 1], (sizeof(double)) * (nX - 1));
                    S->particle_distributions[index + nX - 1] = S->previous_particle_distributions[otherIndex];
                }
            }
        } else if (xDirection == 1) {
            for (int z = 0; z < nZ; z++) {
                int zmd = (nZ + z - zDirection) % nZ;
                int otherZIndex = zmd * nXY + dirIndex;
                int zIndex = z * nXY + dirIndex;
                for (int y = startY; y < endY; y++) {
                    int ymd = y - yDirection;
                    int otherIndex = ymd * nX + otherZIndex;
                    int index = y * nX + zIndex;
                    memcpy(&S->particle_distributions[index + 1], &S->previous_particle_distributions[otherIndex], (sizeof(double)) * (nX - 1));
                    S->particle_distributions[index] = S->previous_particle_distributions[otherIndex + nX - 1];
                }
            }
        }
    }
}




///////////////////////////////////////////////////////////////////
//                   ARRAY VERSIONS                              //             
//////////////////////////////////////////////////////////////////



/*ARRAY BASELINE

Flops: 5 * (nX * nZ * q / 3)
Intops: nX * nZ * q * (39 * nY + 32 / 3)
*/
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
