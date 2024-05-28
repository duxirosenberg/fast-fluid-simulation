#include "LBM.h"
#include <immintrin.h>
/**
 * These functions simulate the collision as described in step 3 of the Time Step Algorithm described in
 * Section 3.3.2 of The Lattice Boltzmann Method book.
 */

// Flops: 26
// Intops: 10
static double calculate_feq_arrays(int index,
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


// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
void collision_arrays(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
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
                    double feq = 2+calculate_feq_arrays(feqIndex, c_s, density_field, velocity_field, (double) directions[3 * i], (double) directions[3 * i + 1], (double) directions[3 * i + 2], weights[i]);

                    int index = x + y * nX + z * nX * nY + i * nX * nY * nZ;
                    //Equation 3.9
                    previous_particle_distributions[index] = omtauinv * particle_distributions[index] + tauinv * feq;

                }
            }
        }
    }
}



// Flops: 26
// Intops: 10
static double calculate_feq_struct( int feqIndex,
                                    struct LBMarrays* S,
                                    int i
)                      {
    double velocityX = S->velocity_field[3 * feqIndex];
    double velocityY = S->velocity_field[3 * feqIndex + 1];
    double velocityZ = S->velocity_field[3 * feqIndex + 2];
    int directionX = S->directions[3 * i];
    int directionY = S->directions[3 * i + 1];
    int directionZ = S->directions[3 * i + 2];
    double weight = S->weights[i];

    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ; // 5 Flops
    // Equation 3.4 with c_s^2 = 1/3
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ; // 5 Flops
    return weight * S->density_field[feqIndex] * (1.0 + dot_product / (S->c_s * S->c_s) + dot_product * dot_product / (2 * S->c_s * S->c_s * S->c_s * S->c_s) - norm_square / (2 * S->c_s * S->c_s)); // 16 Flops
}

// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
void collision_baseline(struct LBMarrays* S) {
    const double tauinv = 1.0 / S->tau;
    const double omtauinv = 1.0 - tauinv;

    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                for(int i = 0; i < S->direction_size; i++) {
                    int feqIndex = (z * S->nX * S->nY) + (y * S->nX) + x;
                    double feq = calculate_feq_struct(feqIndex, S, i); // 10 Intops, 26 Flops

                    int index = x + y * S->nX + z * S->nX * S->nY + i * S->nX * S->nY * S->nZ;
                    // Equation 3.9
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tauinv * feq;
                }
            }
        }
    }
}




// Flops: 21
// Intops: 10
static double calculate_feq_struct2( int feqIndex,
                                    struct LBMarrays* S,
                                    int i
)                      {
    double velocityX = S->velocity_field[3 * feqIndex];
    double velocityY = S->velocity_field[3 * feqIndex + 1];
    double velocityZ = S->velocity_field[3 * feqIndex + 2];
    int directionX = S->directions[3 * i];
    int directionY = S->directions[3 * i + 1];
    int directionZ = S->directions[3 * i + 2];
    double weight = S->weights[i];

    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ; // 5 Flops
    // Equation 3.4 with c_s^2 = 1/3
    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ; // 5 Flops
    return weight * S->density_field[feqIndex] * (1.0 + dot_product / (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2)); // 11 Flops
}



void collision_1(struct LBMarrays* S) {
    const double omtauinv = 1.0 - S->tau_inv;

    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                for(int i = 0; i < S->direction_size; i++) {
                    int feqIndex = (z * S->nXY) + (y * S->nX) + x;
                    double feq = calculate_feq_struct2(feqIndex, S, i); // 10 Intops, 21 Flops

                    int index = x + y * S->nX + z * S->nXY + i * S->nXYZ;
                    // Equation 3.9
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                }
            }
        }
    }
}
// 

// function inlined triple loop reduced to single loope
void collision_2(struct LBMarrays* S) {

        const double omtauinv = 1.0 - S->tau_inv;

        for(int feqIndex = 0; feqIndex < S->nXYZ; ++feqIndex) {

                    double velocityX = S->velocity_field[3 * feqIndex];
                    double velocityY = S->velocity_field[3 * feqIndex + 1];
                    double velocityZ = S->velocity_field[3 * feqIndex + 2];
                    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

            for(int i = 0; i < S->direction_size; i++) {

                //  switched is half as fast for 8 and 32 ..
                    double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));

                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                }
            }   
}


// cache optimised...
// function inlined triple loop reduced to single loope
// need as baseline ...
// cache optimised
// blocking for i; with 3 as a start since 9, 15 and 27 are all divided by 3
// for us new baseline same as above

void collision_3(struct LBMarrays* S) {
        const double omtauinv = 1.0 - S->tau_inv;

            for(int i = 0; i < S->direction_size; i++) {

                //  switched is half as fast for 8 and 32 ..
                    double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = 0; feqIndex < S->nXYZ; ++feqIndex) {

                    double velocityX = S->velocity_field[3 * feqIndex];
                    double velocityY = S->velocity_field[3 * feqIndex + 1];
                    double velocityZ = S->velocity_field[3 * feqIndex + 2];

                    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
                    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));

                    const int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                }
            }
}





void collision_4(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = 256;
    const double omtauinv = 1.0 - S->tau_inv;

        // for(int ii = 0; ii < S->direction_size; ii+=iBB) {
                // for(int i = ii; i < ii+iBB; ++i) {



        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {



                    double velocityX = S->velocity_field[3 * feqIndex];
                    double velocityY = S->velocity_field[3 * feqIndex + 1];
                    double velocityZ = S->velocity_field[3 * feqIndex + 2];


                    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                
            }
        }
    }
}


void collision_5(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = 512;
    const double omtauinv = 1.0 - S->tau_inv;

        // for(int ii = 0; ii < S->direction_size; ii+=iBB) {
                // for(int i = ii; i < ii+iBB; ++i) {



        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {



                    double velocityX = S->velocity_field[3 * feqIndex];
                    double velocityY = S->velocity_field[3 * feqIndex + 1];
                    double velocityZ = S->velocity_field[3 * feqIndex + 2];


                    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                
            }
        }
    }
}

// blocking with bigger ...
void collision_7(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = 2048;
    const double omtauinv = 1.0 - S->tau_inv;

        // for(int ii = 0; ii < S->direction_size; ii+=iBB) {
                // for(int i = ii; i < ii+iBB; ++i) {



        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {



                    double velocityX = S->velocity_field[3 * feqIndex];
                    double velocityY = S->velocity_field[3 * feqIndex + 1];
                    double velocityZ = S->velocity_field[3 * feqIndex + 2];


                    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                
            }
        }
    }
}

void collision_6(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = 512;
    const double omtauinv = 1.0 - S->tau_inv;

        // for(int ii = 0; ii < S->direction_size; ii+=iBB) {
                // for(int i = ii; i < ii+iBB; ++i) {



        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {



                    double velocityX = S->velocity_field[3 * feqIndex];
                    double velocityY = S->velocity_field[3 * feqIndex + 1];
                    double velocityZ = S->velocity_field[3 * feqIndex + 2];


                    double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                
            }
        }
    }
}

// first try without testing error :D:D

void collision_8(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = 1024;
    const double omtauinv = 1.0 - S->tau_inv;

// 0xFFFFFFFF : 0
// __m256i _mm256_cmpeq_epi32 (__m256i a, __m256i b)
// __m128i _mm_cmpeq_epi32 (__m128i a, __m128i b)
// __m128i _mm_set_epi32 (int e3, int e2, int e1, int e0)
// __m256i _mm256_set_epi32 (int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0)


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        // __m128i ones4_i32 =  _mm_set_epi32(1, 1, 1, 1);
        // __m128i true4_i32 = _mm_cmpeq_epi32 (ones4_i32, ones4_i32);


        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {

                    double weight = S->weights[i];



                    // make 4 in a row on 3 vectors...
                    // __m128i dirs_tmp = _mm_loadu_si128(S->directions + 3*i);
                    __m128i dirs_tmp = _mm_maskload_epi32 (S->directions + 3*i, true4_i32);
                    __m256d dirsX  = _mm256_cvtepi32_pd(dirs_tmp);




                    // const int directionX = S->directions[3 * i];
                    // const int directionY = S->directions[3 * i + 1];
                    // const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {



                    __m256d velsX = _mm256_loadu_pd(S->velocity_field + 3*i);
                    // double velocityX = S->velocity_field[3 * feqIndex];
                    // double velocityY = S->velocity_field[3 * feqIndex + 1];
                    // double velocityZ = S->velocity_field[3 * feqIndex + 2];


                    // double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    // double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
                    // _MM_SHUFFLE(3,2,1,0)  not shuffle ..

                    __m256d dv_X2   = _mm256_mul_pd(dirsX, velsX);
                    __m256d velsX2 = _mm256_mul_pd(velsX, velsX);



                    __m256d dv_X_r0    = _mm256_hadd_pd(dv_X2, dv_X2);
                    __m256d vels_X2_r0 = _mm256_hadd_pd(velsX2, velsX2);


                    __m256d dv_X_r1    = _mm256_permute4x64_pd(dv_X2, _MM_SHUFFLE(1,0,3,2));
                    __m256d vels_X2_r1 = _mm256_permute4x64_pd(velsX2, _MM_SHUFFLE(1,0,3,2));

                    __m256d dv_dot     =_mm256_add_pd (dv_X_r0, dv_X_r1);
                    __m256d vels_norm  =_mm256_add_pd (vels_X2_r0, vels_X2_r1);

                    double dot_product = _mm256_cvtsd_f64(dv_dot);
                    double norm_square = _mm256_cvtsd_f64(vels_norm);

                    

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product / 
                    (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                
            }
        }
    }
}






void collision_9(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = 1024;
    const double omtauinv = 1.0 - S->tau_inv;


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);

        __m256d ones_pd = _mm256_set1_pd(1.0);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/S->c_s2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/S->c_s4);
    
        __m256d TAUinv_pd = _mm256_set1_pd(S->tau_inv);
        __m256d omTAUinv_pd = _mm256_set1_pd(omtauinv);
        


        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i+=4) {

                    double weight = S->weights[i];
                    __m256d weight_4 = _mm256_set1_pd (weight);


                    __m128i dirs_tmp0 = _mm_maskload_epi32 (S->directions + 3*i, true4_i32);
                    __m128i dirs_tmp1 = _mm_maskload_epi32 (S->directions + 3*i+4, true4_i32);
                    __m128i dirs_tmp2 = _mm_maskload_epi32 (S->directions + 3*i+8, true4_i32);

                    __m256d dirsX_0  = _mm256_cvtepi32_pd(dirs_tmp0);
                    __m256d dirsX_1  = _mm256_cvtepi32_pd(dirs_tmp1);
                    __m256d dirsX_2  = _mm256_cvtepi32_pd(dirs_tmp2);


                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {



                    __m256d vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*i);
                    __m256d vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*i+ 4);
                    __m256d vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*i+ 8);

                    //////////////////////////////////////////////////////

                    __m256d dv_X2_0   = _mm256_mul_pd(dirsX_0, vels_X_0);
                    __m256d dv_X2_1   = _mm256_mul_pd(dirsX_1, vels_X_1);
                    __m256d dv_X2_2   = _mm256_mul_pd(dirsX_2, vels_X_2);
                    
                    __m256d vel_X2_0 = _mm256_mul_pd(vels_X_0, vels_X_0);
                    __m256d vel_X2_1 = _mm256_mul_pd(vels_X_1, vels_X_1);
                    __m256d vel_X2_2=  _mm256_mul_pd(vels_X_2, vels_X_2);

                    /////////////////////////////////

                    // _MM_SHUFFLE(0,2,0,1) = 0b00100001

                    // due to o3 compiler and ssa i dont think we have to reorder for latencies ...

                    __m256d dv_A    = _mm256_hadd_pd(dv_X2_0, dv_X2_1); //7
                    __m256d dv_B    = _mm256_hadd_pd(dv_X2_1,dv_X2_2);  //7
                    __m256d dv_rest = _mm256_permute2f128_pd(dv_X2_0, dv_X2_2, _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d dv_AB   = _mm256_blend_pd(dv_A,dv_B, 0b1100); //2 
                    __m256d dv_dot  = _mm256_add_pd(dv_AB,dv_rest);  // 4
                    

                    __m256d vel_A    = _mm256_hadd_pd(vel_X2_0, vel_X2_1); //7
                    __m256d vel_B    = _mm256_hadd_pd(vel_X2_1,vel_X2_2); //7
                    __m256d vel_rest = _mm256_permute2f128_pd(vel_X2_0, vel_X2_2, _MM_SHUFFLE(0,2,0,1));//3
                    __m256d vel_AB   = _mm256_blend_pd(vel_A,vel_B, 0b1100); //2
                    __m256d v_norm   = _mm256_add_pd(vel_AB,vel_rest);//4


                    // double dot_product = _mm256_cvtsd_f64(dv_dot);
                    // double norm_square = _mm256_cvtsd_f64(vels_norm);
                    int index = feqIndex + i * S->nXYZ;
                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);


                    // maybe unroll 4 times to avoid latencies, not sure if compiler does it here but if not its bad ...

                    __m256d A              = _mm256_fmadd_pd (dv_dot,twoCS4inv_pd, twoCS2inv_pd);
                    __m256d B              = _mm256_fmadd_pd (dv_dot, A, ones_pd);
                    __m256d C              = _mm256_fmsub_pd (v_norm, twoCS2inv_pd, B); //-1

                    __m256d preF0          = _mm256_mul_pd(weight_4,df_pd);
                    __m256d preF1          = _mm256_mul_pd(TAUinv_pd,preF0);
                    __m256d F              = _mm256_mul_pd(preF1,C);
                    __m256d RET            = _mm256_fmsub_pd(omTAUinv_pd, PD_pd, C); //-1

                    _mm256_store_pd(S->previous_particle_distributions+index, RET);


                    // const double feq = 
                    // weight * S->density_field[feqIndex] * (1.0+dot_product*cs2     +         dot_product * dot_product*cs4        -      norm_square *cs2;
                    // S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
                
            }
        }
    }
}





// cache optimised
// void collision_3b(struct LBMarrays* S) {
//     const double omtauinv = 1.0 - S->tau_inv;

//         for(int i = 0; i < S->direction_size; i++) {
//                     double weight = S->weights[i];
//                     const int directionX = S->directions[3 * i];
//                     const int directionY = S->directions[3 * i + 1];
//                     const int directionZ = S->directions[3 * i + 2];


//                 for(int z = 0; z < S->nZ; z++) {
//                 for(int y = 0; y < S->nY; y++) {
//                 for(int x = 0; x < S->nX; x++) {

//                     int feqIndex = (z * S->nXY) + (y * S->nX) + x;

//                     double velocityX = S->velocity_field[3 * feqIndex];
//                     double velocityY = S->velocity_field[3 * feqIndex + 1];
//                     double velocityZ = S->velocity_field[3 * feqIndex + 2];


//                     double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
//                     double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

//                     const double feq = 
//                     weight * S->density_field[feqIndex] * (1.0 + dot_product / 
//                     (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));

//                     int index = feqIndex + i * S->nXYZ;
//                     S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
//                 }
//             }
//         }
//     }
// }





























// void collision_5...(struct LBMarrays* S) {
//     const int iBB = 9;
//     const double omtauinv = 1.0 - S->tau_inv;

//         for(int ii = 0; ii < S->direction_size; ii+=iBB) {


//             for(int z = 0; z < S->nZ; z++) {
//             for(int y = 0; y < S->nY; y++) {
//             for(int x = 0; x < S->nX; x++) {

//                 for(int i = ii; i < ii+iBB; ++i) {

//                     double weight = S->weights[i];
//                     const int directionX = S->directions[3 * i];
//                     const int directionY = S->directions[3 * i + 1];
//                     const int directionZ = S->directions[3 * i + 2];

//                     int feqIndex = (z * S->nXY) + (y * S->nX) + x;


//                     double velocityX = S->velocity_field[3 * feqIndex];
//                     double velocityY = S->velocity_field[3 * feqIndex + 1];
//                     double velocityZ = S->velocity_field[3 * feqIndex + 2];


//                     double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
//                     double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

//                     const double feq = 
//                     weight * S->density_field[feqIndex] * (1.0 + dot_product / 
//                     (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



//                     int index = feqIndex + i * S->nXYZ;
//                     S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
//                 }
//             }}}
//         }
    
// }


// void collision_6...(struct LBMarrays* S) {
//     const double omtauinv = 1.0 - S->tau_inv;



//                 for(int z = 0; z < S->nZ; z++) {
//                 for(int y = 0; y < S->nY; y++) {
//                 for(int x = 0; x < S->nX; x++) {
//         for(int i = 0; i < S->direction_size; i++) {
//                     double weight = S->weights[i];
//                     const int directionX = S->directions[3 * i];
//                     const int directionY = S->directions[3 * i + 1];
//                     const int directionZ = S->directions[3 * i + 2];

//                     int feqIndex = (z * S->nXY) + (y * S->nX) + x;


//                     double velocityX = S->velocity_field[3 * feqIndex];
//                     double velocityY = S->velocity_field[3 * feqIndex + 1];
//                     double velocityZ = S->velocity_field[3 * feqIndex + 2];


//                     double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
//                     double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

//                     const double feq = 
//                     weight * S->density_field[feqIndex] * (1.0 + dot_product / 
//                     (S->c_s2) + dot_product * dot_product / (2 * S->c_s4) - norm_square / (2 * S->c_s2));



//                     int index = feqIndex + i * S->nXYZ;
//                     S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + S->tau_inv * feq;
//                 }
//             }
//         }
//     }
// }













// void collision_3_fail(struct LBMarrays* S) {
//     const double omtauinv = 1.0 - S->tau_inv;

//         for(int i = 0; i < S->direction_size; i++) {
//                     double weight = S->weights[i];
//                     const int directionX = S->directions[3 * i];
//                     const int directionY = S->directions[3 * i + 1];
//                     const int directionZ = S->directions[3 * i + 2];


//                 for(int z = 0; z < S->nZ; z++) {
//                 for(int y = 0; y < S->nY; y++) {
//                 for(int x = 0; x < S->nX; x+=2) {

//                     int feqIndex0 = (z * S->nXY) + (y * S->nX) + x;
//                     int feqIndex1 = (z * S->nXY) + (y * S->nX) + x + 1;


//                     double velocityX0 = S->velocity_field[3 * feqIndex0];
//                     double velocityY0 = S->velocity_field[3 * feqIndex0 + 1];
//                     double velocityZ0 = S->velocity_field[3 * feqIndex0 + 2];


//                     double velocityX1 = S->velocity_field[3 * feqIndex1];
//                     double velocityY1 = S->velocity_field[3 * feqIndex1 + 1];
//                     double velocityZ1 = S->velocity_field[3 * feqIndex1 + 2];


//                     double dot_product0 = velocityX0 * directionX + velocityY0 * directionY + velocityZ0 * directionZ;
//                     double norm_square0 = velocityX0 * velocityX0 + velocityY0 *  velocityY0 + velocityZ0 * velocityZ0;


//                     double dot_product1 = velocityX1 * directionX + velocityY1 * directionY + velocityZ1 * directionZ;
//                     double norm_square1 = velocityX1 * velocityX1 + velocityY1 * velocityY1 + velocityZ1 * velocityZ1;

//                     const double feq0 = 
//                     weight * S->density_field[feqIndex0] * (1.0 + dot_product0 / 
//                     (S->c_s2) + dot_product0 * dot_product0 / (2 * S->c_s4) - norm_square0 / (2 * S->c_s2));


//                     const double feq1 = 
//                     weight * S->density_field[feqIndex1] * (1.0 + dot_product1 / 
//                     (S->c_s2) + dot_product1 * dot_product1 / (2 * S->c_s4) - norm_square1 / (2 * S->c_s2));

//                     int index0 = feqIndex0 + i * S->nXYZ;
//                     int index1 = feqIndex1 + i * S->nXYZ;
                    
//                     S->previous_particle_distributions[index0] = omtauinv * S->particle_distributions[index0] + S->tau_inv * feq0;
//                     S->previous_particle_distributions[index1] = omtauinv * S->particle_distributions[index1] + S->tau_inv * feq1;
//                 }
//             }
//         }
//     }
// }