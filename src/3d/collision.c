#include "LBM.h"
#include <immintrin.h>


#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdbool.h>

#define BLOCKSIZE_COL 256


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
static double calculate_feq_struct1( int feqIndex,
                                    struct LBMarrays* S,
                                    int i
)                      {

    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
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
    return weight * S->density_field[feqIndex] * (1.0 + dot_product / (cs2) + dot_product * dot_product / (2 * cs4) - norm_square / (2 * cs2)); // 11 Flops
}



void collision_1(struct LBMarrays* S) {
    const double omtauinv = 1.0 - 1/S->tau;

    for(int x = 0; x < S->nX; x++) {
        for(int y = 0; y < S->nY; y++) {
            for(int z = 0; z < S->nZ; z++) {
                for(int i = 0; i < S->direction_size; i++) {
                    int feqIndex = (z * S->nXY) + (y * S->nX) + x;
                    double feq = calculate_feq_struct1(feqIndex, S, i); // 10 Intops, 21 Flops

                    int index = x + y * S->nX + z * S->nXY + i * S->nXYZ;
                    // Equation 3.9
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + 1/S->tau * feq;
                }
            }
        }
    }
}



void collision_2(struct LBMarrays* S) {

        const double cs2 = S->c_s*S->c_s;
        const double cs4 = cs2*cs2;
        const double tau_inv =  1.0/S->tau;
        const double cs2_inv =  1.0 / cs2;
        const double half_cs4_inv = 0.5 / cs4;
        
        const double omtauinv = 1.0 - tau_inv;
        const double half_cs2_inv = 0.5 * cs2_inv;


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
                    weight * S->density_field[feqIndex] * (1.0 + dot_product* cs2_inv
                    + dot_product * dot_product * half_cs4_inv - norm_square * half_cs2_inv);

                    int index = feqIndex + i * S->nXYZ;
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tau_inv * feq;
                }
            }   
}

void collision_3(struct LBMarrays* S) {
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
    const double tau_inv =  1.0/S->tau;
    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;;
    const double half_cs4_inv = 0.5 / cs4;
        
    for(int i = 0; i < S->direction_size; i++) {
        //  switched is half as fast for 8 and 32 ..
        const double weight = S->weights[i];
        const int directionX = S->directions[3 * i];
        const int directionY = S->directions[3 * i + 1];
        const int directionZ = S->directions[3 * i + 2];

        for(int feqIndex = 0; feqIndex < S->nXYZ; ++feqIndex) {

            const int index = feqIndex + i * S->nXYZ;
            const double velocityX = S->velocity_field[3 * feqIndex];
            const double velocityY = S->velocity_field[3 * feqIndex + 1];
            const double velocityZ = S->velocity_field[3 * feqIndex + 2];

            const double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
            const double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product* cs2_inv
                    + dot_product * dot_product * half_cs4_inv - norm_square * half_cs2_inv);

            S->particle_distributions[index] = omtauinv * S->particle_distributions[index] + tau_inv * feq;
        }
    }
    double* tmp = S->previous_particle_distributions;
    S->previous_particle_distributions = S->particle_distributions;
    S->particle_distributions = tmp;
}



void collision_6(struct LBMarrays* S) {
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;

    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
    const double tau_inv =  1.0/S->tau;
    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;;
    const double half_cs4_inv = 0.5 / cs4;

    const int rest = S->nXYZ % nBB;
    const int limit0 = S->nXYZ - rest;

    for(int nxyz = 0; nxyz < limit0; nxyz+=nBB) {
        for(int i = 0; i < S->direction_size; i++) {
            const double weight = S->weights[i];
            const int directionX = S->directions[3 * i];
            const int directionY = S->directions[3 * i + 1];
            const int directionZ = S->directions[3 * i + 2];

            for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {
                const double velocityX = S->velocity_field[3 * feqIndex];
                const double velocityY = S->velocity_field[3 * feqIndex + 1];
                const double velocityZ = S->velocity_field[3 * feqIndex + 2];


                const double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                const double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product* cs2_inv
                    + dot_product * dot_product * half_cs4_inv - norm_square * half_cs2_inv);

                    int index = feqIndex + i * S->nXYZ;
                    S->particle_distributions[index] = omtauinv * S->particle_distributions[index] + tau_inv * feq;
                
            }
        }
    }
        for(int i = 0; i < S->direction_size; i++) {
                    const double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

            for(int feqIndex = limit0; feqIndex <  S->nXYZ; ++feqIndex) {
                    const double velocityX = S->velocity_field[3 * feqIndex];
                    const double velocityY = S->velocity_field[3 * feqIndex + 1];
                    const double velocityZ = S->velocity_field[3 * feqIndex + 2];
                    const double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    const double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;
                    const double feq = 
                    weight * S->density_field[feqIndex] * (1.0 + dot_product* cs2_inv
                    + dot_product * dot_product * half_cs4_inv - norm_square * half_cs2_inv);
                    int index = feqIndex + i * S->nXYZ;
                    S->particle_distributions[index] = omtauinv * S->particle_distributions[index] + tau_inv * feq;
            }
        }  

    double* tmp = S->previous_particle_distributions;
    S->previous_particle_distributions = S->particle_distributions;
    S->particle_distributions = tmp;
}





void collision_SSA3(struct LBMarrays* S) {
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
    const double tau_inv =  1.0/S->tau;
    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;
    const double half_cs4_inv = 0.5 / cs4;

    const int rest = S->nXYZ % nBB;
    const int limit0 = S->nXYZ - rest;


    for(int nxyz = 0; nxyz < limit0; nxyz+=nBB) {
        for(int i = 0; i < S->direction_size; i++) {
            const double weight = S->weights[i];
            const int directionX = S->directions[3 * i];
            const int directionY = S->directions[3 * i + 1];
            const int directionZ = S->directions[3 * i + 2];

            for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {
                const int index = feqIndex + i * S->nXYZ;

                const double velocityX = S->velocity_fieldX[feqIndex];
                const double velocityY = S->velocity_fieldY[feqIndex];
                const double velocityZ = S->velocity_fieldZ[feqIndex];

                const double vdx = velocityX * directionX;
                const double vdy = velocityY * directionY;
                const double vdz = velocityZ * directionZ;
                const double vdxy         = vdx + vdy ;
                const double dot_product  = vdxy + vdz;

                const double vvx = velocityX * velocityX;
                const double vvy = velocityY * velocityY;
                const double vvz = velocityZ * velocityZ;
                const double vvxy         = vvx + vvy;
                const double norm_square  = vvxy + vvz;


                const double A      = weight * S->density_field[feqIndex];
                const double B      = dot_product * dot_product;
                const double C      = norm_square * half_cs2_inv;
                const double D      = dot_product*cs2_inv;
                const double E      = B * half_cs4_inv;
                const double F      = 1.0 + D;
                const double G      = F +   E;
                const double H      = G -   C;
                const double feq    = A *   H;
                const double I      = tau_inv * feq;
                const double J      = omtauinv*S->particle_distributions[index];
                S->particle_distributions[index] = I+J;
            }
        }
    }
        for(int i = 0; i < S->direction_size; i++) {
                    const double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = limit0; feqIndex < S->nXYZ; ++feqIndex) {

                    
                    const int index = feqIndex + i * S->nXYZ;                    

                    const double velocityX = S->velocity_fieldX[feqIndex];
                    const double velocityY = S->velocity_fieldY[feqIndex];
                    const double velocityZ = S->velocity_fieldZ[feqIndex];

                    const double vdx = velocityX * directionX;
                    const double vdy = velocityY * directionY; 
                    const double vdz = velocityZ * directionZ;
                    const double vdxy         = vdx + vdy ;
                    const double dot_product  = vdxy +vdz;

                    const double vvx = velocityX * velocityX;
                    const double vvy = velocityY * velocityY; 
                    const double vvz = velocityZ * velocityZ;
                    const double vvxy         = vvx + vvy;
                    const double norm_square  = vvxy +vvz;


                    const double A      = weight * S->density_field[feqIndex];
                    const double B      = dot_product * dot_product;
                    const double C      = norm_square * half_cs2_inv;
                    const double D      = dot_product*cs2_inv;
                    const double E      = B *   half_cs4_inv;
                    const double F      = 1.0 + D;
                    const double G      = F +   E;
                    const double H      = G -   C;
                    const double feq    = A *   H;
                    const double I      = tau_inv * feq;
                    const double J      = omtauinv*S->particle_distributions[index];
                    S->particle_distributions[index] = I+J;                
            }
        }

    double* tmp = S->previous_particle_distributions;
    S->previous_particle_distributions = S->particle_distributions;
    S->particle_distributions = tmp;
}


void collision_SSA3_nb(struct LBMarrays* S) {
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
    const double tau_inv =  1.0/S->tau;
    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;
    const double half_cs4_inv = 0.5 / cs4;


        for(int i = 0; i < S->direction_size; i++) {
            const double weight = S->weights[i];
            const int directionX = S->directions[3 * i];
            const int directionY = S->directions[3 * i + 1];
            const int directionZ = S->directions[3 * i + 2];

            for(int feqIndex = 0; feqIndex < S->nXYZ; ++feqIndex) {
                const int index = feqIndex + i * S->nXYZ;

                const double velocityX = S->velocity_fieldX[feqIndex];
                const double velocityY = S->velocity_fieldY[feqIndex];
                const double velocityZ = S->velocity_fieldZ[feqIndex];

                const double vdx = velocityX * directionX;
                const double vdy = velocityY * directionY;
                const double vdz = velocityZ * directionZ;
                const double vdxy         = vdx + vdy ;
                const double dot_product  = vdxy + vdz;

                const double vvx = velocityX * velocityX;
                const double vvy = velocityY * velocityY;
                const double vvz = velocityZ * velocityZ;
                const double vvxy         = vvx + vvy;
                const double norm_square  = vvxy + vvz;


                const double A      = weight * S->density_field[feqIndex];
                const double B      = dot_product * dot_product;
                const double C      = norm_square * half_cs2_inv;
                const double D      = dot_product*cs2_inv;
                const double E      = B * half_cs4_inv;
                const double F      = 1.0 + D;
                const double G      = F +   E;
                const double H      = G -   C;
                const double feq    = A *   H;
                const double I      = tau_inv * feq;
                const double J      = omtauinv*S->particle_distributions[index];
                S->particle_distributions[index] = I+J;
            }
        }
    

    double* tmp = S->previous_particle_distributions;
    S->previous_particle_distributions = S->particle_distributions;
    S->particle_distributions = tmp;
}


void collision_AVX4(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d twos_pd      = _mm256_set1_pd(2.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d halfCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d halfCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);

        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    double weight = S->weights[i];

                    __m256d weight_4 = _mm256_set1_pd(weight);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);

                    __m256d dirX = _mm256_set1_pd(S->directions[3 * i]);
                    __m256d dirY = _mm256_set1_pd(S->directions[3 * i+1]);
                    __m256d dirZ = _mm256_set1_pd(S->directions[3 * i+2]);

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=4) {

                    
                    const int index = feqIndex + i * S->nXYZ;
                    const int b_index = feqIndex + i * S->nXYZ;                    
                    
                    

                    __m256d velX = _mm256_loadu_pd(S->velocity_fieldX + feqIndex);
                    __m256d velY = _mm256_loadu_pd(S->velocity_fieldY + feqIndex);
                    __m256d velZ = _mm256_loadu_pd(S->velocity_fieldZ + feqIndex);
                                                                                     


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);
                    

                    __m256d dotX     = _mm256_mul_pd(   velX, dirX);
                    __m256d dotXY    = _mm256_fmadd_pd( velY, dirY, dotX);
                    __m256d dv_dot    = _mm256_fmadd_pd(velZ, dirZ, dotXY);

                    __m256d normX     = _mm256_mul_pd(  velX, velX);
                    __m256d normXY    = _mm256_fmadd_pd(velY, velY, normX);
                    __m256d v_norm    = _mm256_fmadd_pd(velZ, velZ, normXY);


                    __m256d tmpB          = _mm256_mul_pd(df_pd, tmpA);
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd, PD_pd);

                    // __m256d A              = _mm256_fmadd_pd(dv_dot,twoCS2inv_pd, ones_pd);
                    // __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    // __m256d C              = _mm256_fmadd_pd(B, twoCS2inv_pd, ones_pd);
                    // __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);
                    __m256d A              = _mm256_fmadd_pd(dv_dot,halfCS2inv_pd, twos_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    __m256d C              = _mm256_fmadd_pd(B, halfCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);
                    

                
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
            }
        }
    }
}



void collision_AVX5_u2(struct LBMarrays* S){

    const int nBB = S->nXYZ < BLOCKSIZE_COL ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;

    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;
    const double half_cs4_inv = 0.5 / cs4;

    const int rest = S->nXYZ % BLOCKSIZE_COL;
    const int limit0 = S->nXYZ - rest;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);
        

        for(int nxyz = 0; nxyz < limit0; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    __m256d weight_4 = _mm256_set1_pd(S->weights[i]);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);

                    __m256d dirX = _mm256_set1_pd(S->directions[3 * i]);
                    __m256d dirY = _mm256_set1_pd(S->directions[3 * i+1]);
                    __m256d dirZ = _mm256_set1_pd(S->directions[3 * i+2]);

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=16) {

                    
                    const int b_feqIndex = feqIndex + 4;
                    const int c_feqIndex = feqIndex + 8;
                    const int d_feqIndex = feqIndex + 12;
                    const int index      = feqIndex   + i * S->nXYZ;
                    const int b_index    = b_feqIndex + i * S->nXYZ;    
                    const int c_index    = c_feqIndex + i * S->nXYZ;    
                    const int d_index    = d_feqIndex + i * S->nXYZ;      

                    // LOADING      
                    
                    __m256d velX = _mm256_loadu_pd(S->velocity_fieldX + feqIndex);
                    __m256d velY = _mm256_loadu_pd(S->velocity_fieldY + feqIndex);
                    __m256d velZ = _mm256_loadu_pd(S->velocity_fieldZ + feqIndex);

                    __m256d b_velX = _mm256_loadu_pd(S->velocity_fieldX + b_feqIndex);
                    __m256d b_velY = _mm256_loadu_pd(S->velocity_fieldY + b_feqIndex);
                    __m256d b_velZ = _mm256_loadu_pd(S->velocity_fieldZ + b_feqIndex);

                    __m256d c_velX = _mm256_loadu_pd(S->velocity_fieldX + c_feqIndex);
                    __m256d c_velY = _mm256_loadu_pd(S->velocity_fieldY + c_feqIndex);
                    __m256d c_velZ = _mm256_loadu_pd(S->velocity_fieldZ + c_feqIndex);

                    __m256d d_velX = _mm256_loadu_pd(S->velocity_fieldX + d_feqIndex);
                    __m256d d_velY = _mm256_loadu_pd(S->velocity_fieldY + d_feqIndex);
                    __m256d d_velZ = _mm256_loadu_pd(S->velocity_fieldZ + d_feqIndex);




                    __m256d df_pd = _mm256_load_pd(S->density_field          + feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions + index);

                    __m256d b_df_pd = _mm256_load_pd(S->density_field          + b_feqIndex);
                    __m256d b_PD_pd = _mm256_load_pd(S->particle_distributions + b_index);

                    __m256d c_df_pd = _mm256_load_pd(S->density_field          + c_feqIndex);
                    __m256d c_PD_pd = _mm256_load_pd(S->particle_distributions + c_index);

                    __m256d d_df_pd = _mm256_load_pd(S->density_field          + d_feqIndex);
                    __m256d d_PD_pd = _mm256_load_pd(S->particle_distributions + d_index);


                    // COMPUTATION DOTPROUCT & NORM

                    __m256d dotX     = _mm256_mul_pd(   velX, dirX);
                    __m256d dotXY    = _mm256_fmadd_pd( velY, dirY, dotX);
                    __m256d dv_dot    = _mm256_fmadd_pd(velZ, dirZ, dotXY);
                    __m256d normX     = _mm256_mul_pd(  velX, velX);
                    __m256d normXY    = _mm256_fmadd_pd(velY, velY, normX);
                    __m256d v_norm    = _mm256_fmadd_pd(velZ, velZ, normXY);
                    __m256d b_dotX     = _mm256_mul_pd(   b_velX, dirX);

                    __m256d b_dotXY    = _mm256_fmadd_pd( b_velY, dirY,   b_dotX);
                    __m256d b_dv_dot    = _mm256_fmadd_pd(b_velZ, dirZ,   b_dotXY);
                    __m256d b_normX     = _mm256_mul_pd(  b_velX, b_velX);
                    __m256d b_normXY    = _mm256_fmadd_pd(b_velY, b_velY, b_normX);
                    __m256d b_v_norm    = _mm256_fmadd_pd(b_velZ, b_velZ, b_normXY);

                    __m256d c_dotX     = _mm256_mul_pd(   c_velX, dirX);
                    __m256d c_dotXY    = _mm256_fmadd_pd( c_velY, dirY,   c_dotX);
                    __m256d c_dv_dot    = _mm256_fmadd_pd(c_velZ, dirZ,   c_dotXY);
                    __m256d c_normX     = _mm256_mul_pd(  c_velX, c_velX);
                    __m256d c_normXY    = _mm256_fmadd_pd(c_velY, c_velY, c_normX);
                    __m256d c_v_norm    = _mm256_fmadd_pd(c_velZ, c_velZ, c_normXY);

                    __m256d d_dotX     = _mm256_mul_pd(   d_velX, dirX);
                    __m256d d_dotXY    = _mm256_fmadd_pd( d_velY, dirY,   d_dotX);
                    __m256d d_dv_dot    = _mm256_fmadd_pd(d_velZ, dirZ,   d_dotXY);
                    __m256d d_normX     = _mm256_mul_pd(  d_velX, d_velX);
                    __m256d d_normXY    = _mm256_fmadd_pd(d_velY, d_velY, d_normX);
                    __m256d d_v_norm    = _mm256_fmadd_pd(d_velZ, d_velZ, d_normXY);


                    // //////////////////////////////
                    // COMPUTATION RETURN
                    // //////////////////////////////

                    __m256d tmpB          = _mm256_mul_pd(    df_pd,       tmpA);
                    __m256d tmpC          = _mm256_mul_pd(    omTAUinv_pd, PD_pd);
                    __m256d dotSqr        = _mm256_mul_pd(dv_dot, dv_dot);
                    __m256d A             = _mm256_fmsub_pd( dotSqr,      CS2inv_pd, v_norm);
                    __m256d B             = _mm256_fmadd_pd(dv_dot, CS2inv_pd, ones_pd);
                    __m256d C             = _mm256_fmadd_pd(A, twoCS2inv_pd, B);
                    __m256d D             = _mm256_fmadd_pd( tmpB,        C,            tmpC);


                    __m256d b_tmpB          = _mm256_mul_pd(    b_df_pd,       tmpA);
                    __m256d b_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   b_PD_pd);
                    __m256d b_dotSqr        = _mm256_mul_pd(b_dv_dot, b_dv_dot);
                    __m256d b_A             = _mm256_fmsub_pd( b_dotSqr,      CS2inv_pd, b_v_norm);
                    __m256d b_B             = _mm256_fmadd_pd(b_dv_dot, CS2inv_pd, ones_pd);
                    __m256d b_C             = _mm256_fmadd_pd(b_A, twoCS2inv_pd, b_B);
                    __m256d b_D             = _mm256_fmadd_pd( b_tmpB,        b_C,            b_tmpC);


                    __m256d c_tmpB          = _mm256_mul_pd(    c_df_pd,       tmpA);
                    __m256d c_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   c_PD_pd);
                    __m256d c_dotSqr        = _mm256_mul_pd(c_dv_dot, c_dv_dot);
                    __m256d c_A             = _mm256_fmsub_pd( c_dotSqr,      CS2inv_pd, c_v_norm);
                    __m256d c_B             = _mm256_fmadd_pd(c_dv_dot, CS2inv_pd, ones_pd);
                    __m256d c_C             = _mm256_fmadd_pd(c_A, twoCS2inv_pd, c_B);
                    __m256d c_D             = _mm256_fmadd_pd( c_tmpB,        c_C,            c_tmpC);


                    __m256d d_tmpB          = _mm256_mul_pd(    d_df_pd,       tmpA);
                    __m256d d_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   d_PD_pd);
                    __m256d d_dotSqr        = _mm256_mul_pd(d_dv_dot, d_dv_dot);
                    __m256d d_A             = _mm256_fmsub_pd( d_dotSqr,      CS2inv_pd, d_v_norm);
                    __m256d d_B             = _mm256_fmadd_pd(d_dv_dot, CS2inv_pd, ones_pd);
                    __m256d d_C             = _mm256_fmadd_pd(d_A, twoCS2inv_pd, d_B);
                    __m256d d_D             = _mm256_fmadd_pd( d_tmpB,        d_C,            d_tmpC);           

                    //////////////////
                    // STORE
                    //////////////////
                    
                    
                    _mm256_storeu_pd(S->particle_distributions+index,       D);
                    _mm256_storeu_pd(S->particle_distributions+b_index,   b_D);
                    _mm256_storeu_pd(S->particle_distributions+c_index,   c_D);
                    _mm256_storeu_pd(S->particle_distributions+d_index,   d_D);
            }
        }
    }
        for(int i = 0; i < S->direction_size; i++) {
                    const double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = limit0; feqIndex < S->nXYZ; ++feqIndex) {

                    
                    const int index = feqIndex + i * S->nXYZ;                    

                    const double velocityX = S->velocity_fieldX[feqIndex];
                    const double velocityY = S->velocity_fieldY[feqIndex];
                    const double velocityZ = S->velocity_fieldZ[feqIndex];

                    const double vdx = velocityX * directionX;
                    const double vdy = velocityY * directionY; 
                    const double vdz = velocityZ * directionZ;
                    const double vdxy         = vdx + vdy ;
                    const double dot_product  = vdxy +vdz;

                    const double vvx = velocityX * velocityX;
                    const double vvy = velocityY * velocityY; 
                    const double vvz = velocityZ * velocityZ;
                    const double vvxy         = vvx + vvy;
                    const double norm_square  = vvxy +vvz;


                    const double A      = weight * S->density_field[feqIndex];
                    const double B      = dot_product * dot_product;
                    const double C      = norm_square * half_cs2_inv;
                    const double D      = dot_product*cs2_inv;
                    const double E      = B *   half_cs4_inv;
                    const double F      = 1.0 + D;
                    const double G      = F +   E;
                    const double H      = G -   C;
                    const double feq    = A *   H;
                    const double I      = tau_inv * feq;
                    const double J      = omtauinv*S->particle_distributions[index];
                    S->particle_distributions[index] = I+J;                
            }
        }


    double* tmp = S->previous_particle_distributions;
    S->previous_particle_distributions = S->particle_distributions;
    S->particle_distributions = tmp;
}



void collision_AVX5_u2_nb(struct LBMarrays* S){

    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;

    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;
    const double half_cs4_inv = 0.5 / cs4;
    const int rest = S->nXYZ % 16;
    const int limit0 = S->nXYZ - rest;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);
        

            for(int i = 0; i < S->direction_size; i++) { 

                    __m256d weight_4 = _mm256_set1_pd(S->weights[i]);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);

                    __m256d dirX = _mm256_set1_pd(S->directions[3 * i]);
                    __m256d dirY = _mm256_set1_pd(S->directions[3 * i+1]);
                    __m256d dirZ = _mm256_set1_pd(S->directions[3 * i+2]);

                for(int feqIndex = 0; feqIndex < limit0; feqIndex+=16) {

                    
                    const int b_feqIndex = feqIndex + 4;
                    const int c_feqIndex = feqIndex + 8;
                    const int d_feqIndex = feqIndex + 12;
                    const int index      = feqIndex   + i * S->nXYZ;
                    const int b_index    = b_feqIndex + i * S->nXYZ;    
                    const int c_index    = c_feqIndex + i * S->nXYZ;    
                    const int d_index    = d_feqIndex + i * S->nXYZ;      

                    // LOADING      
                    
                    __m256d velX = _mm256_loadu_pd(S->velocity_fieldX + feqIndex);
                    __m256d velY = _mm256_loadu_pd(S->velocity_fieldY + feqIndex);
                    __m256d velZ = _mm256_loadu_pd(S->velocity_fieldZ + feqIndex);

                    __m256d b_velX = _mm256_loadu_pd(S->velocity_fieldX + b_feqIndex);
                    __m256d b_velY = _mm256_loadu_pd(S->velocity_fieldY + b_feqIndex);
                    __m256d b_velZ = _mm256_loadu_pd(S->velocity_fieldZ + b_feqIndex);

                    __m256d c_velX = _mm256_loadu_pd(S->velocity_fieldX + c_feqIndex);
                    __m256d c_velY = _mm256_loadu_pd(S->velocity_fieldY + c_feqIndex);
                    __m256d c_velZ = _mm256_loadu_pd(S->velocity_fieldZ + c_feqIndex);

                    __m256d d_velX = _mm256_loadu_pd(S->velocity_fieldX + d_feqIndex);
                    __m256d d_velY = _mm256_loadu_pd(S->velocity_fieldY + d_feqIndex);
                    __m256d d_velZ = _mm256_loadu_pd(S->velocity_fieldZ + d_feqIndex);




                    __m256d df_pd = _mm256_load_pd(S->density_field          + feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions + index);

                    __m256d b_df_pd = _mm256_load_pd(S->density_field          + b_feqIndex);
                    __m256d b_PD_pd = _mm256_load_pd(S->particle_distributions + b_index);

                    __m256d c_df_pd = _mm256_load_pd(S->density_field          + c_feqIndex);
                    __m256d c_PD_pd = _mm256_load_pd(S->particle_distributions + c_index);

                    __m256d d_df_pd = _mm256_load_pd(S->density_field          + d_feqIndex);
                    __m256d d_PD_pd = _mm256_load_pd(S->particle_distributions + d_index);


                    // COMPUTATION DOTPROUCT & NORM

                    __m256d dotX     = _mm256_mul_pd(   velX, dirX);
                    __m256d dotXY    = _mm256_fmadd_pd( velY, dirY, dotX);
                    __m256d dv_dot    = _mm256_fmadd_pd(velZ, dirZ, dotXY);
                    __m256d normX     = _mm256_mul_pd(  velX, velX);
                    __m256d normXY    = _mm256_fmadd_pd(velY, velY, normX);
                    __m256d v_norm    = _mm256_fmadd_pd(velZ, velZ, normXY);
                    __m256d b_dotX     = _mm256_mul_pd(   b_velX, dirX);

                    __m256d b_dotXY    = _mm256_fmadd_pd( b_velY, dirY,   b_dotX);
                    __m256d b_dv_dot    = _mm256_fmadd_pd(b_velZ, dirZ,   b_dotXY);
                    __m256d b_normX     = _mm256_mul_pd(  b_velX, b_velX);
                    __m256d b_normXY    = _mm256_fmadd_pd(b_velY, b_velY, b_normX);
                    __m256d b_v_norm    = _mm256_fmadd_pd(b_velZ, b_velZ, b_normXY);

                    __m256d c_dotX     = _mm256_mul_pd(   c_velX, dirX);
                    __m256d c_dotXY    = _mm256_fmadd_pd( c_velY, dirY,   c_dotX);
                    __m256d c_dv_dot    = _mm256_fmadd_pd(c_velZ, dirZ,   c_dotXY);
                    __m256d c_normX     = _mm256_mul_pd(  c_velX, c_velX);
                    __m256d c_normXY    = _mm256_fmadd_pd(c_velY, c_velY, c_normX);
                    __m256d c_v_norm    = _mm256_fmadd_pd(c_velZ, c_velZ, c_normXY);

                    __m256d d_dotX     = _mm256_mul_pd(   d_velX, dirX);
                    __m256d d_dotXY    = _mm256_fmadd_pd( d_velY, dirY,   d_dotX);
                    __m256d d_dv_dot    = _mm256_fmadd_pd(d_velZ, dirZ,   d_dotXY);
                    __m256d d_normX     = _mm256_mul_pd(  d_velX, d_velX);
                    __m256d d_normXY    = _mm256_fmadd_pd(d_velY, d_velY, d_normX);
                    __m256d d_v_norm    = _mm256_fmadd_pd(d_velZ, d_velZ, d_normXY);


                    // //////////////////////////////
                    // COMPUTATION RETURN
                    // //////////////////////////////

                    __m256d tmpB          = _mm256_mul_pd(    df_pd,       tmpA);
                    __m256d tmpC          = _mm256_mul_pd(    omTAUinv_pd, PD_pd);
                    __m256d dotSqr        = _mm256_mul_pd(dv_dot, dv_dot);
                    __m256d A             = _mm256_fmsub_pd( dotSqr,      CS2inv_pd, v_norm);
                    __m256d B             = _mm256_fmadd_pd(dv_dot, CS2inv_pd, ones_pd);
                    __m256d C             = _mm256_fmadd_pd(A, twoCS2inv_pd, B);
                    __m256d D             = _mm256_fmadd_pd( tmpB,        C,            tmpC);


                    __m256d b_tmpB          = _mm256_mul_pd(    b_df_pd,       tmpA);
                    __m256d b_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   b_PD_pd);
                    __m256d b_dotSqr        = _mm256_mul_pd(b_dv_dot, b_dv_dot);
                    __m256d b_A             = _mm256_fmsub_pd( b_dotSqr,      CS2inv_pd, b_v_norm);
                    __m256d b_B             = _mm256_fmadd_pd(b_dv_dot, CS2inv_pd, ones_pd);
                    __m256d b_C             = _mm256_fmadd_pd(b_A, twoCS2inv_pd, b_B);
                    __m256d b_D             = _mm256_fmadd_pd( b_tmpB,        b_C,            b_tmpC);


                    __m256d c_tmpB          = _mm256_mul_pd(    c_df_pd,       tmpA);
                    __m256d c_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   c_PD_pd);
                    __m256d c_dotSqr        = _mm256_mul_pd(c_dv_dot, c_dv_dot);
                    __m256d c_A             = _mm256_fmsub_pd( c_dotSqr,      CS2inv_pd, c_v_norm);
                    __m256d c_B             = _mm256_fmadd_pd(c_dv_dot, CS2inv_pd, ones_pd);
                    __m256d c_C             = _mm256_fmadd_pd(c_A, twoCS2inv_pd, c_B);
                    __m256d c_D             = _mm256_fmadd_pd( c_tmpB,        c_C,            c_tmpC);


                    __m256d d_tmpB          = _mm256_mul_pd(    d_df_pd,       tmpA);
                    __m256d d_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   d_PD_pd);
                    __m256d d_dotSqr        = _mm256_mul_pd(d_dv_dot, d_dv_dot);
                    __m256d d_A             = _mm256_fmsub_pd( d_dotSqr,      CS2inv_pd, d_v_norm);
                    __m256d d_B             = _mm256_fmadd_pd(d_dv_dot, CS2inv_pd, ones_pd);
                    __m256d d_C             = _mm256_fmadd_pd(d_A, twoCS2inv_pd, d_B);
                    __m256d d_D             = _mm256_fmadd_pd( d_tmpB,        d_C,            d_tmpC);           

                    //////////////////
                    // STORE
                    //////////////////
                    
                    
                    _mm256_storeu_pd(S->particle_distributions+index,       D);
                    _mm256_storeu_pd(S->particle_distributions+b_index,   b_D);
                    _mm256_storeu_pd(S->particle_distributions+c_index,   c_D);
                    _mm256_storeu_pd(S->particle_distributions+d_index,   d_D);
            }

            const double weight = S->weights[i];
            const int directionX = S->directions[3 * i];
            const int directionY = S->directions[3 * i + 1];
            const int directionZ = S->directions[3 * i + 2];

            for(int feqIndex = limit0; feqIndex < S->nXYZ; ++feqIndex) {

                    
                    const int index = feqIndex + i * S->nXYZ;                    

                    const double velocityX = S->velocity_fieldX[feqIndex];
                    const double velocityY = S->velocity_fieldY[feqIndex];
                    const double velocityZ = S->velocity_fieldZ[feqIndex];

                    const double vdx = velocityX * directionX;
                    const double vdy = velocityY * directionY; 
                    const double vdz = velocityZ * directionZ;
                    const double vdxy         = vdx + vdy ;
                    const double dot_product  = vdxy +vdz;

                    const double vvx = velocityX * velocityX;
                    const double vvy = velocityY * velocityY; 
                    const double vvz = velocityZ * velocityZ;
                    const double vvxy         = vvx + vvy;
                    const double norm_square  = vvxy +vvz;


                    const double A      = weight * S->density_field[feqIndex];
                    const double B      = dot_product * dot_product;
                    const double C      = norm_square * half_cs2_inv;
                    const double D      = dot_product*cs2_inv;
                    const double E      = B *   half_cs4_inv;
                    const double F      = 1.0 + D;
                    const double G      = F +   E;
                    const double H      = G -   C;
                    const double feq    = A *   H;
                    const double I      = tau_inv * feq;
                    const double J      = omtauinv*S->particle_distributions[index];
                    S->particle_distributions[index] = I+J;                
            }
            
        }
    


    double* tmp = S->previous_particle_distributions;
    S->previous_particle_distributions = S->particle_distributions;
    S->particle_distributions = tmp;
}



// BLOCKING TEST ...
// choice between 256, 512, 1024 seems to be not relevant


        // __m256d twos_pd      = _mm256_set1_pd(2.0);

                    // __m256d A              = _mm256_fmadd_pd( dv_dot,      twoCS2inv_pd, twos_pd);
                    // __m256d B              = _mm256_fmsub_pd( dv_dot,      A,            v_norm);
                    // __m256d C              = _mm256_fmadd_pd( B,           twoCS2inv_pd, ones_pd);
                    // __m256d D              = _mm256_fmadd_pd( tmpB,        C,            tmpC);
                    // __m256d b_A              = _mm256_fmadd_pd( b_dv_dot,      twoCS2inv_pd,   twos_pd);
                    // __m256d b_B              = _mm256_fmsub_pd( b_dv_dot,      b_A,            b_v_norm);
                    // __m256d b_C              = _mm256_fmadd_pd( b_B,           twoCS2inv_pd,   ones_pd);
                    // __m256d b_D              = _mm256_fmadd_pd( b_tmpB,        b_C,            b_tmpC);
                    // __m256d c_A              = _mm256_fmadd_pd( c_dv_dot,      twoCS2inv_pd,   ones_pd);
                    // __m256d c_B              = _mm256_fmsub_pd( c_dv_dot,      c_A,            c_v_norm);
                    // __m256d c_C              = _mm256_fmadd_pd( c_B,           twoCS2inv_pd,   ones_pd);
                    // __m256d c_D              = _mm256_fmadd_pd( c_tmpB,        c_C,            c_tmpC);
                    // __m256d d_A              = _mm256_fmadd_pd( d_dv_dot,      twoCS2inv_pd,   twos_pd);
                    // __m256d d_B              = _mm256_fmsub_pd( d_dv_dot,      d_A,            d_v_norm);
                    // __m256d d_C              = _mm256_fmadd_pd( d_B,           twoCS2inv_pd,   ones_pd);
                    // __m256d d_D              = _mm256_fmadd_pd( d_tmpB,        d_C,            d_tmpC);


