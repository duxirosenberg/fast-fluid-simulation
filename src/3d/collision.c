#include "LBM.h"
#include <immintrin.h>


#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdbool.h>

#define BLOCKSIZE_COL 256
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

            S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tau_inv * feq;
        }
    }
}



void collision_6(struct LBMarrays* S) {
    // const int iBB = 9;
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;

    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
    const double tau_inv =  1.0/S->tau;
    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;;
    const double half_cs4_inv = 0.5 / cs4;

    for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
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
                    S->previous_particle_distributions[index] = omtauinv * S->particle_distributions[index] + tau_inv * feq;
                
            }
        }
    }
}



void collision_SSA(struct LBMarrays* S) {


        const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
        const double cs2 = S->c_s*S->c_s;
        const double cs4 = cs2*cs2;
        const double tau_inv =  1.0/S->tau;
        const double omtauinv = 1.0 - tau_inv;
        const double cs2_inv =  1.0 / cs2;
        const double half_cs2_inv = 0.5 * cs2_inv;;
        const double half_cs4_inv = 0.5 / cs4;


        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    const double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; ++feqIndex) {

                    
                    const int index = feqIndex + i * S->nXYZ;
                    // const int b_index = feqIndex + i * S->nXYZ;
                    

                    const double velocityX = S->velocity_field[3 * feqIndex];
                    const double velocityY = S->velocity_field[3 * feqIndex + 1];
                    const double velocityZ = S->velocity_field[3 * feqIndex + 2];


                    const double dot_product = velocityX * directionX + velocityY * directionY + velocityZ * directionZ;
                    const double norm_square = velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ;

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
                    S->previous_particle_distributions[index] = I+J;

                
            }
        }
    }
}



void collision_SSA2(struct LBMarrays* S) {
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;
    const double tau_inv =  1.0/S->tau;
    const double omtauinv = 1.0 - tau_inv;
    const double cs2_inv =  1.0 / cs2;
    const double half_cs2_inv = 0.5 * cs2_inv;
    const double half_cs4_inv = 0.5 / cs4;


    for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
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
                S->previous_particle_distributions[index] = I+J;
            }
        }
    }
}





void collision_SSA_u(struct LBMarrays* S) {
    

        const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
        const double cs2 = S->c_s*S->c_s;
        const double cs4 = cs2*cs2;
        const double tau_inv =  1.0/S->tau;
        const double omtauinv = 1.0 - tau_inv;
        const double cs2_inv =  1.0 / cs2;
        const double half_cs2_inv = 0.5 * cs2_inv;;
        const double half_cs4_inv = 0.5 / cs4;


        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    const double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int a_feqIndex = nxyz; a_feqIndex < nxyz+nBB; a_feqIndex+=2) {


                    
                    const int b_feqIndex = a_feqIndex + 1;
                    const int a_index = a_feqIndex + i * S->nXYZ;
                    const int b_index = b_feqIndex + i * S->nXYZ;



                    const double a_velocityX = S->velocity_fieldX[a_feqIndex];
                    const double a_velocityY = S->velocity_fieldY[a_feqIndex];
                    const double a_velocityZ = S->velocity_fieldZ[a_feqIndex];

                    const double b_velocityX = S->velocity_fieldX[b_feqIndex];
                    const double b_velocityY = S->velocity_fieldY[b_feqIndex];
                    const double b_velocityZ = S->velocity_fieldZ[b_feqIndex];



                    const double a_dot_product = a_velocityX * directionX + a_velocityY * directionY + a_velocityZ * directionZ;
                    const double a_norm_square = a_velocityX * a_velocityX + a_velocityY * a_velocityY + a_velocityZ * a_velocityZ;
                    const double b_dot_product = b_velocityX * directionX + b_velocityY * directionY + b_velocityZ * directionZ;
                    const double b_norm_square = b_velocityX * b_velocityX + b_velocityY * b_velocityY + b_velocityZ * b_velocityZ;

                    const double a_A      = weight * S->density_field[a_feqIndex];
                    const double b_A      = weight * S->density_field[b_feqIndex];
                    const double a_B      = a_dot_product * a_dot_product;
                    const double b_B      = b_dot_product * b_dot_product;
                    const double a_C      = a_norm_square * half_cs2_inv;
                    const double b_C      = b_norm_square * half_cs2_inv;
                    const double a_D      = a_dot_product*  cs2_inv;
                    const double b_D      = b_dot_product*  cs2_inv;
                    const double a_E      = a_B *           half_cs4_inv;
                    const double b_E      = b_B *           half_cs4_inv;
                    const double a_F      = 1.0 + a_D;
                    const double b_F      = 1.0 + b_D;
                    const double a_G      = a_F +   a_E;
                    const double b_G      = b_F +   b_E;
                    const double a_H      = a_G -   a_C;
                    const double b_H      = b_G -   b_C;
                    const double a_feq    = a_A *   a_H;
                    const double b_feq    = b_A *   b_H;
                    const double a_I      = tau_inv * a_feq;
                    const double b_I      = tau_inv * b_feq;
                    const double a_J      = omtauinv*S->particle_distributions[a_index];
                    const double b_J      = omtauinv*S->particle_distributions[b_index];
                    const double a_ppd    = a_I+a_J;
                    const double b_ppd    = b_I+b_J;
                    S->previous_particle_distributions[a_index] = a_ppd;
                    S->previous_particle_distributions[b_index] = b_ppd;

                
            }
        }
    }
}



void collision_SSA2_u(struct LBMarrays* S) {

        const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
        const double cs2 = S->c_s*S->c_s;
        const double cs4 = cs2*cs2;
        const double tau_inv =  1.0/S->tau;
        const double omtauinv = 1.0 - tau_inv;
        const double cs2_inv =  1.0 / cs2;
        const double half_cs2_inv = 0.5 * cs2_inv;;
        const double half_cs4_inv = 0.5 / cs4;


        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) {
                    const double weight = S->weights[i];
                    const int directionX = S->directions[3 * i];
                    const int directionY = S->directions[3 * i + 1];
                    const int directionZ = S->directions[3 * i + 2];

                for(int a_feqIndex = nxyz; a_feqIndex < nxyz+nBB; a_feqIndex+=2) {
                    
                    const int b_feqIndex = a_feqIndex + 1;
                    const int a_index = a_feqIndex + i * S->nXYZ;
                    const int b_index = b_feqIndex + i * S->nXYZ;

                    const double a_velocityX = S->velocity_field[3 * a_feqIndex];
                    const double a_velocityY = S->velocity_field[3 * a_feqIndex + 1];
                    const double a_velocityZ = S->velocity_field[3 * a_feqIndex + 2];
                    const double b_velocityX = S->velocity_field[3 * b_feqIndex];
                    const double b_velocityY = S->velocity_field[3 * b_feqIndex + 1];
                    const double b_velocityZ = S->velocity_field[3 * b_feqIndex + 2];

                    const double a_dot_product = a_velocityX * directionX + a_velocityY * directionY + a_velocityZ * directionZ;
                    const double a_norm_square = a_velocityX * a_velocityX + a_velocityY * a_velocityY + a_velocityZ * a_velocityZ;
                    const double b_dot_product = b_velocityX * directionX + b_velocityY * directionY + b_velocityZ * directionZ;
                    const double b_norm_square = b_velocityX * b_velocityX + b_velocityY * b_velocityY + b_velocityZ * b_velocityZ;

                    const double a_A      = weight * S->density_field[a_feqIndex];
                    const double a_B      = a_dot_product * a_dot_product;
                    const double a_C      = a_norm_square * half_cs2_inv;
                    const double a_D      = a_dot_product*  cs2_inv;
                    const double a_E      = a_B *           half_cs4_inv;
                    const double a_F      = 1.0 + a_D;
                    const double a_G      = a_F +   a_E;
                    const double a_H      = a_G -   a_C;
                    const double a_feq    = a_A *   a_H;
                    const double a_I      = tau_inv * a_feq;
                    const double a_J      = omtauinv*S->particle_distributions[a_index];
                    const double a_ppd    = a_I+a_J;

                    const double b_A      = weight * S->density_field[b_feqIndex];
                    const double b_B      = b_dot_product * b_dot_product;
                    const double b_C      = b_norm_square * half_cs2_inv;
                    const double b_D      = b_dot_product*  cs2_inv;
                    const double b_E      = b_B *           half_cs4_inv;
                    const double b_F      = 1.0 + b_D;
                    const double b_G      = b_F +   b_E;
                    const double b_H      = b_G -   b_C;
                    const double b_feq    = b_A *   b_H;
                    const double b_I      = tau_inv * b_feq;
                    const double b_J      = omtauinv*S->particle_distributions[b_index];
                    const double b_ppd    = b_I+b_J;
                    S->previous_particle_distributions[a_index] = a_ppd;
                    S->previous_particle_distributions[b_index] = b_ppd;

                
            }
        }
    }
}



void collision_AVX(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    double weight = S->weights[i];

                    __m256d weight_4 = _mm256_set1_pd(weight);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);


                    __m128i dirs_tmp = _mm_maskload_epi32(S->directions + 3*i, true4_i32);
                    __m256d dirsX  = _mm256_cvtepi32_pd(dirs_tmp);
                    __m256d dirsX_0 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(0,2,1,0));
                    __m256d dirsX_1 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(1,0,2,1));
                    __m256d dirsX_2 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(2,1,0,2));

        
            // for(int feqIndex = 0; feqIndex < S->nXYZ; feqIndex+=4) {
            for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=4) {


                    int index = feqIndex + i * S->nXYZ;

                    __m256d vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex);
                    __m256d vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 4);
                    __m256d vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 8);


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);


                    __m256d dv_X2_0   = _mm256_mul_pd(dirsX_0, vels_X_0);
                    __m256d dv_X2_1   = _mm256_mul_pd(dirsX_1, vels_X_1);
                    __m256d dv_X2_2   = _mm256_mul_pd(dirsX_2, vels_X_2);
                    
                    __m256d vel_X2_0 = _mm256_mul_pd(vels_X_0, vels_X_0);
                    __m256d vel_X2_1 = _mm256_mul_pd(vels_X_1, vels_X_1);
                    __m256d vel_X2_2=  _mm256_mul_pd(vels_X_2, vels_X_2);

                    // //////////////////////////////////////////////

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


                    __m256d tmpB          = _mm256_mul_pd(df_pd, tmpA);
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd, PD_pd);

                    __m256d A              = _mm256_fmadd_pd(dv_dot,twoCS2inv_pd, ones_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    __m256d C              = _mm256_fmadd_pd(B, twoCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);

                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
            }
        }
    }
}



void collision_AVX_u(struct LBMarrays* S) {
    
        const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;

        const double cs2 = S->c_s*S->c_s;
        const double cs4 = cs2*cs2;
        const double tau_inv = 1/S->tau;
        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    double weight = S->weights[i];

                    __m256d weight_4 = _mm256_set1_pd(weight);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);


                    __m128i dirs_tmp = _mm_maskload_epi32(S->directions + 3*i, true4_i32);
                    __m256d dirsX  = _mm256_cvtepi32_pd(dirs_tmp);
                    __m256d dirsX_0 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(0,2,1,0));
                    __m256d dirsX_1 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(1,0,2,1));
                    __m256d dirsX_2 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(2,1,0,2));

            // unroll once
            for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=8) {


                    int index = feqIndex + i * S->nXYZ;
                    int b_feqIndex = feqIndex + 4;
                    int b_index = index + 4;

                    // a ....

                    __m256d vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex);
                    __m256d vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 4);
                    __m256d vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 8);


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);


                    __m256d dv_X2_0   = _mm256_mul_pd(dirsX_0, vels_X_0);
                    __m256d dv_X2_1   = _mm256_mul_pd(dirsX_1, vels_X_1);
                    __m256d dv_X2_2   = _mm256_mul_pd(dirsX_2, vels_X_2);
                    
                    __m256d vel_X2_0 = _mm256_mul_pd(vels_X_0, vels_X_0);
                    __m256d vel_X2_1 = _mm256_mul_pd(vels_X_1, vels_X_1);
                    __m256d vel_X2_2 = _mm256_mul_pd(vels_X_2, vels_X_2);

                    // b ...


                    __m256d b_vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*b_feqIndex);
                    __m256d b_vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*b_feqIndex+ 4);
                    __m256d b_vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*b_feqIndex+ 8);


                    __m256d b_df_pd = _mm256_load_pd(S->density_field+b_feqIndex);
                    __m256d b_PD_pd = _mm256_load_pd(S->particle_distributions+b_index);


                    __m256d b_dv_X2_0   = _mm256_mul_pd(dirsX_0, b_vels_X_0);
                    __m256d b_dv_X2_1   = _mm256_mul_pd(dirsX_1, b_vels_X_1);
                    __m256d b_dv_X2_2   = _mm256_mul_pd(dirsX_2, b_vels_X_2);
                    
                    __m256d b_vel_X2_0 = _mm256_mul_pd(b_vels_X_0, b_vels_X_0);
                    __m256d b_vel_X2_1 = _mm256_mul_pd(b_vels_X_1, b_vels_X_1);
                    __m256d b_vel_X2_2 = _mm256_mul_pd(b_vels_X_2, b_vels_X_2);



                    // //////////////////////////////////////////////
                    // //////////////////////////////////////////////

                    // a ....

                    __m256d dv_A     = _mm256_hadd_pd( dv_X2_0,  dv_X2_1); //7
                    __m256d dv_B     = _mm256_hadd_pd( dv_X2_1,  dv_X2_2);  //7
                    __m256d dv_AB    = _mm256_blend_pd(dv_A,     dv_B, 0b1100); //2 
                    __m256d vel_A    = _mm256_hadd_pd( vel_X2_0, vel_X2_1); //7
                    __m256d vel_B    = _mm256_hadd_pd( vel_X2_1, vel_X2_2); //7
                    __m256d vel_AB   = _mm256_blend_pd(vel_A,    vel_B, 0b1100); //2

                    __m256d dv_rest  = _mm256_permute2f128_pd(dv_X2_0,   dv_X2_2,   _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d vel_rest = _mm256_permute2f128_pd(vel_X2_0, vel_X2_2,  _MM_SHUFFLE(0,2,0,1));//3

                    __m256d dv_dot   = _mm256_add_pd(dv_AB,      dv_rest);  // 4
                    __m256d v_norm   = _mm256_add_pd(vel_AB,    vel_rest);//4


                    __m256d tmpB          = _mm256_mul_pd(df_pd,        tmpA);
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd,  PD_pd);

                    __m256d A              = _mm256_fmadd_pd(dv_dot, twoCS2inv_pd,  ones_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A,             v_norm);
                    __m256d C              = _mm256_fmadd_pd(B,      twoCS2inv_pd,  ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,   C,             tmpC);

                    //  b ...

                    __m256d b_dv_A     = _mm256_hadd_pd( b_dv_X2_0,  b_dv_X2_1); //7
                    __m256d b_dv_B     = _mm256_hadd_pd( b_dv_X2_1,  b_dv_X2_2);  //7
                    __m256d b_dv_AB    = _mm256_blend_pd(b_dv_A,     b_dv_B, 0b1100); //2 
                    __m256d b_vel_A    = _mm256_hadd_pd( b_vel_X2_0, b_vel_X2_1); //7
                    __m256d b_vel_B    = _mm256_hadd_pd( b_vel_X2_1, b_vel_X2_2); //7
                    __m256d b_vel_AB   = _mm256_blend_pd(b_vel_A,    b_vel_B, 0b1100); //2

                    __m256d b_dv_rest  = _mm256_permute2f128_pd(b_dv_X2_0,  b_dv_X2_2,   _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d b_vel_rest = _mm256_permute2f128_pd(b_vel_X2_0, b_vel_X2_2,  _MM_SHUFFLE(0,2,0,1));//3

                    __m256d b_dv_dot   = _mm256_add_pd(b_dv_AB,     b_dv_rest);  // 4
                    __m256d b_v_norm   = _mm256_add_pd(b_vel_AB,    b_vel_rest);//4


                    __m256d b_tmpB          = _mm256_mul_pd(b_df_pd,        tmpA);
                    __m256d b_tmpC          = _mm256_mul_pd(omTAUinv_pd,  b_PD_pd);

                    __m256d b_A              = _mm256_fmadd_pd(b_dv_dot, twoCS2inv_pd,    ones_pd);
                    __m256d b_B              = _mm256_fmsub_pd(b_dv_dot, b_A,             b_v_norm);
                    __m256d b_C              = _mm256_fmadd_pd(b_B,      twoCS2inv_pd,    ones_pd);
                    __m256d b_D              = _mm256_fmadd_pd(b_tmpB,   b_C,             b_tmpC);




                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
                    _mm256_storeu_pd(S->previous_particle_distributions+b_index, b_D);
            }
        }
    }
}



void collision_AVX_u2(struct LBMarrays* S) {
    


    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    double weight = S->weights[i];

                    __m256d weight_4 = _mm256_set1_pd(weight);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);


                    __m128i dirs_tmp = _mm_maskload_epi32(S->directions + 3*i, true4_i32);
                    __m256d dirsX  = _mm256_cvtepi32_pd(dirs_tmp);
                    __m256d dirsX_0 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(0,2,1,0));
                    __m256d dirsX_1 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(1,0,2,1));
                    __m256d dirsX_2 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(2,1,0,2));

            // unroll twice
            for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=16) {


                    int index = feqIndex + i * S->nXYZ;
                    int b_feqIndex = feqIndex + 4;
                    int b_index = index +       4;

                    int c_feqIndex = feqIndex + 8;
                    int c_index = index +       8;

                    int d_feqIndex = feqIndex + 12;
                    int d_index = index +       12;

                    // a ....
                    __m256d vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex);
                    __m256d vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 4);
                    __m256d vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 8);


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);


                    __m256d dv_X2_0   = _mm256_mul_pd(dirsX_0, vels_X_0);
                    __m256d dv_X2_1   = _mm256_mul_pd(dirsX_1, vels_X_1);
                    __m256d dv_X2_2   = _mm256_mul_pd(dirsX_2, vels_X_2);
                    
                    __m256d vel_X2_0 = _mm256_mul_pd(vels_X_0, vels_X_0);
                    __m256d vel_X2_1 = _mm256_mul_pd(vels_X_1, vels_X_1);
                    __m256d vel_X2_2 = _mm256_mul_pd(vels_X_2, vels_X_2);

                    // b ...
                    __m256d b_vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*b_feqIndex);
                    __m256d b_vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*b_feqIndex+ 4);
                    __m256d b_vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*b_feqIndex+ 8);


                    __m256d b_df_pd = _mm256_load_pd(S->density_field+b_feqIndex);
                    __m256d b_PD_pd = _mm256_load_pd(S->particle_distributions+b_index);


                    __m256d b_dv_X2_0   = _mm256_mul_pd(dirsX_0, b_vels_X_0);
                    __m256d b_dv_X2_1   = _mm256_mul_pd(dirsX_1, b_vels_X_1);
                    __m256d b_dv_X2_2   = _mm256_mul_pd(dirsX_2, b_vels_X_2);
                    
                    __m256d b_vel_X2_0 = _mm256_mul_pd(b_vels_X_0, b_vels_X_0);
                    __m256d b_vel_X2_1 = _mm256_mul_pd(b_vels_X_1, b_vels_X_1);
                    __m256d b_vel_X2_2 = _mm256_mul_pd(b_vels_X_2, b_vels_X_2);

                    // c ...
                    __m256d c_vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*c_feqIndex);
                    __m256d c_vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*c_feqIndex+ 4);
                    __m256d c_vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*c_feqIndex+ 8);


                    __m256d c_df_pd = _mm256_load_pd(S->density_field+c_feqIndex);
                    __m256d c_PD_pd = _mm256_load_pd(S->particle_distributions+c_index);


                    __m256d c_dv_X2_0   = _mm256_mul_pd(dirsX_0, c_vels_X_0);
                    __m256d c_dv_X2_1   = _mm256_mul_pd(dirsX_1, c_vels_X_1);
                    __m256d c_dv_X2_2   = _mm256_mul_pd(dirsX_2, c_vels_X_2);
                    
                    __m256d c_vel_X2_0 = _mm256_mul_pd(c_vels_X_0, c_vels_X_0);
                    __m256d c_vel_X2_1 = _mm256_mul_pd(c_vels_X_1, c_vels_X_1);
                    __m256d c_vel_X2_2 = _mm256_mul_pd(c_vels_X_2, c_vels_X_2);

                    // d ...
                    __m256d d_vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*d_feqIndex);
                    __m256d d_vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*d_feqIndex+ 4);
                    __m256d d_vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*d_feqIndex+ 8);


                    __m256d d_df_pd = _mm256_load_pd(S->density_field+d_feqIndex);
                    __m256d d_PD_pd = _mm256_load_pd(S->particle_distributions+d_index);


                    __m256d d_dv_X2_0   = _mm256_mul_pd(dirsX_0, d_vels_X_0);
                    __m256d d_dv_X2_1   = _mm256_mul_pd(dirsX_1, d_vels_X_1);
                    __m256d d_dv_X2_2   = _mm256_mul_pd(dirsX_2, d_vels_X_2);
                    
                    __m256d d_vel_X2_0 = _mm256_mul_pd(d_vels_X_0, d_vels_X_0);
                    __m256d d_vel_X2_1 = _mm256_mul_pd(d_vels_X_1, d_vels_X_1);
                    __m256d d_vel_X2_2 = _mm256_mul_pd(d_vels_X_2, d_vels_X_2);



                    // //////////////////////////////////////////////
                    // //////////////////////////////////////////////

                    // a ....
                    __m256d dv_A     = _mm256_hadd_pd( dv_X2_0,  dv_X2_1); //7
                    __m256d dv_B     = _mm256_hadd_pd( dv_X2_1,  dv_X2_2);  //7
                    __m256d dv_AB    = _mm256_blend_pd(dv_A,     dv_B, 0b1100); //2 
                    __m256d vel_A    = _mm256_hadd_pd( vel_X2_0, vel_X2_1); //7
                    __m256d vel_B    = _mm256_hadd_pd( vel_X2_1, vel_X2_2); //7
                    __m256d vel_AB   = _mm256_blend_pd(vel_A,    vel_B, 0b1100); //2

                    __m256d dv_rest  = _mm256_permute2f128_pd(dv_X2_0,   dv_X2_2,   _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d vel_rest = _mm256_permute2f128_pd(vel_X2_0, vel_X2_2,  _MM_SHUFFLE(0,2,0,1));//3

                    __m256d dv_dot   = _mm256_add_pd(dv_AB,      dv_rest);  // 4
                    __m256d v_norm   = _mm256_add_pd(vel_AB,    vel_rest);//4


                    __m256d tmpB          = _mm256_mul_pd(df_pd,        tmpA);//4
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd,  PD_pd);//4

                    __m256d A              = _mm256_fmadd_pd(dv_dot, twoCS2inv_pd,  ones_pd);//4
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A,             v_norm);//4
                    __m256d C              = _mm256_fmadd_pd(B,      twoCS2inv_pd,  ones_pd);//4
                    __m256d D              = _mm256_fmadd_pd(tmpB,   C,             tmpC);//4

                    //  b ...
                    __m256d b_dv_A     = _mm256_hadd_pd( b_dv_X2_0,  b_dv_X2_1); //7
                    __m256d b_dv_B     = _mm256_hadd_pd( b_dv_X2_1,  b_dv_X2_2);  //7
                    __m256d b_dv_AB    = _mm256_blend_pd(b_dv_A,     b_dv_B, 0b1100); //2 
                    __m256d b_vel_A    = _mm256_hadd_pd( b_vel_X2_0, b_vel_X2_1); //7
                    __m256d b_vel_B    = _mm256_hadd_pd( b_vel_X2_1, b_vel_X2_2); //7
                    __m256d b_vel_AB   = _mm256_blend_pd(b_vel_A,    b_vel_B, 0b1100); //2

                    __m256d b_dv_rest  = _mm256_permute2f128_pd(b_dv_X2_0,  b_dv_X2_2,   _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d b_vel_rest = _mm256_permute2f128_pd(b_vel_X2_0, b_vel_X2_2,  _MM_SHUFFLE(0,2,0,1));//3

                    __m256d b_dv_dot   = _mm256_add_pd(b_dv_AB,     b_dv_rest);  // 4
                    __m256d b_v_norm   = _mm256_add_pd(b_vel_AB,    b_vel_rest);//4


                    __m256d b_tmpB          = _mm256_mul_pd(b_df_pd,        tmpA);//4
                    __m256d b_tmpC          = _mm256_mul_pd(omTAUinv_pd,  b_PD_pd);//4

                    __m256d b_A              = _mm256_fmadd_pd(b_dv_dot, twoCS2inv_pd,    ones_pd);//4
                    __m256d b_B              = _mm256_fmsub_pd(b_dv_dot, b_A,             b_v_norm);//4
                    __m256d b_C              = _mm256_fmadd_pd(b_B,      twoCS2inv_pd,    ones_pd);//4
                    __m256d b_D              = _mm256_fmadd_pd(b_tmpB,   b_C,             b_tmpC);//4

                    //c...
                    __m256d c_dv_A     = _mm256_hadd_pd( c_dv_X2_0,  c_dv_X2_1); //7
                    __m256d c_dv_B     = _mm256_hadd_pd( c_dv_X2_1,  c_dv_X2_2);  //7
                    __m256d c_dv_AB    = _mm256_blend_pd(c_dv_A,     c_dv_B, 0b1100); //2 
                    __m256d c_vel_A    = _mm256_hadd_pd( c_vel_X2_0, c_vel_X2_1); //7
                    __m256d c_vel_B    = _mm256_hadd_pd( c_vel_X2_1, c_vel_X2_2); //7
                    __m256d c_vel_AB   = _mm256_blend_pd(c_vel_A,    c_vel_B, 0b1100); //2

                    __m256d c_dv_rest  = _mm256_permute2f128_pd(c_dv_X2_0,  c_dv_X2_2,   _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d c_vel_rest = _mm256_permute2f128_pd(c_vel_X2_0, c_vel_X2_2,  _MM_SHUFFLE(0,2,0,1));//3

                    __m256d c_dv_dot   = _mm256_add_pd(c_dv_AB,     c_dv_rest);  // 4
                    __m256d c_v_norm   = _mm256_add_pd(c_vel_AB,    c_vel_rest);//4


                    __m256d c_tmpB          = _mm256_mul_pd(c_df_pd,        tmpA);//4
                    __m256d c_tmpC          = _mm256_mul_pd(omTAUinv_pd,  c_PD_pd);//4

                    __m256d c_A              = _mm256_fmadd_pd(c_dv_dot, twoCS2inv_pd,    ones_pd);//4
                    __m256d c_B              = _mm256_fmsub_pd(c_dv_dot, c_A,             c_v_norm);//4
                    __m256d c_C              = _mm256_fmadd_pd(c_B,      twoCS2inv_pd,    ones_pd);//4
                    __m256d c_D              = _mm256_fmadd_pd(c_tmpB,   c_C,             c_tmpC);//4

                    //d...
                    __m256d d_dv_A     = _mm256_hadd_pd( d_dv_X2_0,  d_dv_X2_1); //7
                    __m256d d_dv_B     = _mm256_hadd_pd( d_dv_X2_1,  d_dv_X2_2);  //7
                    __m256d d_dv_AB    = _mm256_blend_pd(d_dv_A,     d_dv_B, 0b1100); //2 
                    __m256d d_vel_A    = _mm256_hadd_pd( d_vel_X2_0, d_vel_X2_1); //7
                    __m256d d_vel_B    = _mm256_hadd_pd( d_vel_X2_1, d_vel_X2_2); //7
                    __m256d d_vel_AB   = _mm256_blend_pd(d_vel_A,    d_vel_B, 0b1100); //2

                    __m256d d_dv_rest  = _mm256_permute2f128_pd(d_dv_X2_0,  d_dv_X2_2,   _MM_SHUFFLE(0,2,0,1)); //3
                    __m256d d_vel_rest = _mm256_permute2f128_pd(d_vel_X2_0, d_vel_X2_2,  _MM_SHUFFLE(0,2,0,1));//3

                    __m256d d_dv_dot   = _mm256_add_pd(d_dv_AB,     d_dv_rest);  // 4
                    __m256d d_v_norm   = _mm256_add_pd(d_vel_AB,    d_vel_rest);//4


                    __m256d d_tmpB          = _mm256_mul_pd(d_df_pd,        tmpA);//4
                    __m256d d_tmpC          = _mm256_mul_pd(omTAUinv_pd,  d_PD_pd);//4

                    __m256d d_A              = _mm256_fmadd_pd(d_dv_dot, twoCS2inv_pd,    ones_pd);//4
                    __m256d d_B              = _mm256_fmsub_pd(d_dv_dot, d_A,             d_v_norm);//4
                    __m256d d_C              = _mm256_fmadd_pd(d_B,      twoCS2inv_pd,    ones_pd);//4
                    __m256d d_D              = _mm256_fmadd_pd(d_tmpB,   d_C,             d_tmpC);//4

                    ////////////////////////
                    ////////////////////////


                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index,   D);
                    _mm256_storeu_pd(S->previous_particle_distributions+b_index, b_D);
                    _mm256_storeu_pd(S->previous_particle_distributions+c_index, c_D);
                    _mm256_storeu_pd(S->previous_particle_distributions+d_index, d_D);
            }
        }
    }
}



void collision_AVX2(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    double weight = S->weights[i];

                    __m256d weight_4 = _mm256_set1_pd(weight);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);


                    __m128i dirs_tmp = _mm_maskload_epi32(S->directions + 3*i, true4_i32);
                    __m256d dirsX  = _mm256_cvtepi32_pd(dirs_tmp);
                    __m256d dirsX_0 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(0,2,1,0));
                    __m256d dirsX_1 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(1,0,2,1));
                    __m256d dirsX_2 = _mm256_permute4x64_pd(dirsX, _MM_SHUFFLE(2,1,0,2));

        
            // for(int feqIndex = 0; feqIndex < S->nXYZ; feqIndex+=4) {
            for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=4) {


                    int index = feqIndex + i * S->nXYZ;

                    __m256d vels_X_0 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex);
                    __m256d vels_X_1 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 4);
                    __m256d vels_X_2 = _mm256_loadu_pd(S->velocity_field + 3*feqIndex+ 8);


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);


                    __m256d dv_X2_0   = _mm256_mul_pd(dirsX_0, vels_X_0);
                    __m256d dv_X2_1   = _mm256_mul_pd(dirsX_1, vels_X_1);
                    __m256d dv_X2_2   = _mm256_mul_pd(dirsX_2, vels_X_2);
                    
                    __m256d vel_X2_0 = _mm256_mul_pd(vels_X_0, vels_X_0);
                    __m256d vel_X2_1 = _mm256_mul_pd(vels_X_1, vels_X_1);
                    __m256d vel_X2_2=  _mm256_mul_pd(vels_X_2, vels_X_2);

                    // //////////////////////////////////////////////
                    // try without hadd ...

                    __m256d dv_A0    = _mm256_permute_pd(  dv_X2_0, _MM_SHUFFLE(1,1,1,1)); //1
                    __m256d dv_BC0   = _mm256_permute_pd(  dv_X2_1, _MM_SHUFFLE(1,1,1,1)); //1
                    __m256d dv_D0    = _mm256_permute_pd(  dv_X2_2, _MM_SHUFFLE(1,1,1,1)); //1
                    __m256d dv_AD1   = _mm256_blend_pd(    dv_X2_0, dv_X2_2, 0b1100); //2
                    __m256d dv_AD0   = _mm256_blend_pd(    dv_A0,   dv_D0, 0b1100); //2

                    __m256d vel_A0    = _mm256_permute_pd(  vel_X2_0, _MM_SHUFFLE(1,1,1,1)); //1
                    __m256d vel_BC0   = _mm256_permute_pd(  vel_X2_1, _MM_SHUFFLE(1,1,1,1)); //1
                    __m256d vel_D0    = _mm256_permute_pd(  vel_X2_2, _MM_SHUFFLE(1,1,1,1)); //1
                    __m256d vel_AD1   = _mm256_blend_pd(    vel_X2_0, vel_X2_2, 0b1100); //2
                    __m256d vel_AD0   = _mm256_blend_pd(    vel_A0,   vel_D0, 0b1100); //2



                    __m256d dv_AD    = _mm256_add_pd(      dv_AD0,  dv_AD1);//4
                    __m256d dv_BC    = _mm256_add_pd(      dv_X2_1, dv_BC0);//4
                    __m256d vel_AD    = _mm256_add_pd(      vel_AD0,  vel_AD1);//4
                    __m256d vel_BC    = _mm256_add_pd(      vel_X2_1, vel_BC0);//4



                    __m256d dv_AB0   = _mm256_permute_pd(  dv_X2_0, _MM_SHUFFLE(3,2,3,2)); //3
                    __m256d dv_CD3   = _mm256_permute_pd(  dv_X2_2, _MM_SHUFFLE(1,0,1,0)); //3
                    __m256d dv_ABCD0 = _mm256_blend_pd(    dv_AD,   dv_BC, 0b0110); //2

                    __m256d vel_AB0   = _mm256_permute_pd(  vel_X2_0, _MM_SHUFFLE(3,2,3,2)); //3
                    __m256d vel_CD3   = _mm256_permute_pd(  vel_X2_2, _MM_SHUFFLE(1,0,1,0)); //3
                    __m256d vel_ABCD0 = _mm256_blend_pd(    vel_AD,   vel_BC, 0b0110); //2


                    __m256d dv_ABCD1 = _mm256_blend_pd(    dv_AB0,  dv_CD3, 0b1100); //2
                    __m256d vel_ABCD1 = _mm256_blend_pd(    vel_AB0,  vel_CD3, 0b1100); //2

                    __m256d dv_dot   = _mm256_add_pd(      dv_ABCD0,dv_ABCD1);//4
                    __m256d v_norm   = _mm256_add_pd(      vel_ABCD0,vel_ABCD1);//4
                    


                    __m256d tmpB          = _mm256_mul_pd(df_pd, tmpA);
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd, PD_pd);

                    __m256d A              = _mm256_fmadd_pd(dv_dot,twoCS2inv_pd, ones_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    __m256d C              = _mm256_fmadd_pd(B, twoCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);

                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
            }
        }
    }
}



void collision_AVX3(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        

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
                    

                    __m256d velX = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex],
                                                S->velocity_field[3 * feqIndex+3], 
                                                S->velocity_field[3 * feqIndex+6], 
                                                S->velocity_field[3 * feqIndex+9]
                                                );
                    

                    __m256d velY = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+1],
                                                S->velocity_field[3 * feqIndex+4], 
                                                S->velocity_field[3 * feqIndex+7], 
                                                S->velocity_field[3 * feqIndex+10]
                                                );

                    __m256d velZ = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+2],
                                                S->velocity_field[3 * feqIndex+5], 
                                                S->velocity_field[3 * feqIndex+8], 
                                                S->velocity_field[3 * feqIndex+11]
                                                );
                                                                   


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);
                    

                    __m256d dotX     = _mm256_mul_pd(  velX, dirX);
                    __m256d dotXY    = _mm256_fmadd_pd(velY, dirY, dotX);
                    __m256d dv_dot    = _mm256_fmadd_pd(velZ, dirZ, dotXY);

                    __m256d normX     = _mm256_mul_pd(  velX, velX);
                    __m256d normXY    = _mm256_fmadd_pd(velY, velY, normX);
                    __m256d v_norm    = _mm256_fmadd_pd(velZ, velZ, normXY);



// //////////////////////////////////////////////////////////////////////////////////////////////

                    __m256d tmpB          = _mm256_mul_pd(df_pd, tmpA);
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd, PD_pd);

                    __m256d A              = _mm256_fmadd_pd(dv_dot,twoCS2inv_pd, ones_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    __m256d C              = _mm256_fmadd_pd(B, twoCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);

                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
            }
        }
    }
}



void collision_AVX3_u(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);


        __m128i true4_i32 =  _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    double weight = S->weights[i];

                    __m256d weight_4 = _mm256_set1_pd(weight);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);

                    __m256d dirX = _mm256_set1_pd(S->directions[3 * i]);
                    __m256d dirY = _mm256_set1_pd(S->directions[3 * i+1]);
                    __m256d dirZ = _mm256_set1_pd(S->directions[3 * i+2]);

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=8) {

                    
                    const int index = feqIndex + i * S->nXYZ;
                    
                    const int b_index = feqIndex + 4 + i * S->nXYZ;                    
                    

                    __m256d velX = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex],
                                                S->velocity_field[3 * feqIndex+3], 
                                                S->velocity_field[3 * feqIndex+6], 
                                                S->velocity_field[3 * feqIndex+9]
                                                );
                    

                    __m256d velY = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+1],
                                                S->velocity_field[3 * feqIndex+4], 
                                                S->velocity_field[3 * feqIndex+7], 
                                                S->velocity_field[3 * feqIndex+10]
                                                );

                    __m256d velZ = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+2],
                                                S->velocity_field[3 * feqIndex+5], 
                                                S->velocity_field[3 * feqIndex+8], 
                                                S->velocity_field[3 * feqIndex+11]
                                                );
                                                                   
                    __m256d b_velX = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+12],
                                                S->velocity_field[3 * feqIndex+15], 
                                                S->velocity_field[3 * feqIndex+18], 
                                                S->velocity_field[3 * feqIndex+21]
                                                );
                    

                    __m256d b_velY = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+13],
                                                S->velocity_field[3 * feqIndex+16], 
                                                S->velocity_field[3 * feqIndex+19], 
                                                S->velocity_field[3 * feqIndex+22]
                                                );

                    __m256d b_velZ = _mm256_set_pd(   
                                                S->velocity_field[3 * feqIndex+14],
                                                S->velocity_field[3 * feqIndex+17], 
                                                S->velocity_field[3 * feqIndex+20], 
                                                S->velocity_field[3 * feqIndex+23]
                                                );


                    __m256d df_pd = _mm256_load_pd(S->density_field+feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions+index);

                    __m256d b_df_pd = _mm256_load_pd(S->density_field+feqIndex+4);
                    __m256d b_PD_pd = _mm256_load_pd(S->particle_distributions+b_index);
                    



                    __m256d dotX     = _mm256_mul_pd(   velX, dirX);
                    __m256d dotXY    = _mm256_fmadd_pd( velY, dirY, dotX);
                    __m256d dv_dot    = _mm256_fmadd_pd(velZ, dirZ, dotXY);
                    __m256d normX     = _mm256_mul_pd(  velX, velX);
                    __m256d normXY    = _mm256_fmadd_pd(velY, velY, normX);
                    __m256d v_norm    = _mm256_fmadd_pd(velZ, velZ, normXY);


                    __m256d b_dotX     = _mm256_mul_pd(   b_velX, dirX);
                    __m256d b_dotXY    = _mm256_fmadd_pd( b_velY, dirY, b_dotX);
                    __m256d b_dv_dot    = _mm256_fmadd_pd(b_velZ, dirZ, b_dotXY);
                    __m256d b_normX     = _mm256_mul_pd(  b_velX, b_velX);
                    __m256d b_normXY    = _mm256_fmadd_pd(b_velY, b_velY, b_normX);
                    __m256d b_v_norm    = _mm256_fmadd_pd(b_velZ, b_velZ, b_normXY);



// //////////////////////////////////////////////////////////////////////////////////////////////

                    __m256d tmpB          = _mm256_mul_pd(df_pd, tmpA);
                    __m256d tmpC          = _mm256_mul_pd(omTAUinv_pd, PD_pd);

                    __m256d A              = _mm256_fmadd_pd(dv_dot,twoCS2inv_pd, ones_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    __m256d C              = _mm256_fmadd_pd(B, twoCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);


                    __m256d b_tmpB          = _mm256_mul_pd(b_df_pd,        tmpA);
                    __m256d b_tmpC          = _mm256_mul_pd(omTAUinv_pd,  b_PD_pd);

                    __m256d b_A              = _mm256_fmadd_pd(b_dv_dot, twoCS2inv_pd,    ones_pd);
                    __m256d b_B              = _mm256_fmsub_pd(b_dv_dot, b_A,             b_v_norm);
                    __m256d b_C              = _mm256_fmadd_pd(b_B,      twoCS2inv_pd,    ones_pd);
                    __m256d b_D              = _mm256_fmadd_pd(b_tmpB,   b_C,             b_tmpC);

                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
                    _mm256_storeu_pd(S->previous_particle_distributions+b_index, b_D);
            }
        }
    }
}


void collision_AVX4(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
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

                    __m256d A              = _mm256_fmadd_pd(dv_dot,twoCS2inv_pd, ones_pd);
                    __m256d B              = _mm256_fmsub_pd(dv_dot, A, v_norm);
                    __m256d C              = _mm256_fmadd_pd(B, twoCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd(tmpB,C,tmpC);

                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
            }
        }
    }
}



void collision_AVX4_u(struct LBMarrays* S) {

    // const int nBB = 1024;512
    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);
       

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
            for(int i = 0; i < S->direction_size; i++) { 

                    __m256d weight_4 = _mm256_set1_pd(S->weights[i]);
                    __m256d tmpA     = _mm256_mul_pd(weight_4,TAUinv_pd);

                    __m256d dirX = _mm256_set1_pd(S->directions[3 * i]);
                    __m256d dirY = _mm256_set1_pd(S->directions[3 * i+1]);
                    __m256d dirZ = _mm256_set1_pd(S->directions[3 * i+2]);

                for(int feqIndex = nxyz; feqIndex < nxyz+nBB; feqIndex+=8) {

                    
                    const int b_feqIndex = feqIndex + 4;
                    const int index      = feqIndex   + i * S->nXYZ;
                    const int b_index    = b_feqIndex + i * S->nXYZ;                    
                    
                    

                    __m256d velX = _mm256_loadu_pd(S->velocity_fieldX + feqIndex);
                    __m256d velY = _mm256_loadu_pd(S->velocity_fieldY + feqIndex);
                    __m256d velZ = _mm256_loadu_pd(S->velocity_fieldZ + feqIndex);

                    __m256d b_velX = _mm256_loadu_pd(S->velocity_fieldX + b_feqIndex);
                    __m256d b_velY = _mm256_loadu_pd(S->velocity_fieldY + b_feqIndex);
                    __m256d b_velZ = _mm256_loadu_pd(S->velocity_fieldZ + b_feqIndex);


                                                                                     


                    __m256d df_pd = _mm256_load_pd(S->density_field          + feqIndex);
                    __m256d PD_pd = _mm256_load_pd(S->particle_distributions + index);

                    __m256d b_df_pd = _mm256_load_pd(S->density_field          + b_feqIndex);
                    __m256d b_PD_pd = _mm256_load_pd(S->particle_distributions + b_index);

                    

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


                
                    // ///////////////////////////
                    

                    __m256d tmpB          = _mm256_mul_pd(    df_pd,       tmpA);
                    __m256d tmpC          = _mm256_mul_pd(    omTAUinv_pd, PD_pd);
                    __m256d A              = _mm256_fmadd_pd( dv_dot,      twoCS2inv_pd, ones_pd);
                    __m256d B              = _mm256_fmsub_pd( dv_dot,      A,            v_norm);
                    __m256d C              = _mm256_fmadd_pd( B,           twoCS2inv_pd, ones_pd);
                    __m256d D              = _mm256_fmadd_pd( tmpB,        C,            tmpC);


                    __m256d b_tmpB          = _mm256_mul_pd(    b_df_pd,       tmpA);
                    __m256d b_tmpC          = _mm256_mul_pd(    omTAUinv_pd,   b_PD_pd);
                    __m256d b_A              = _mm256_fmadd_pd( b_dv_dot,      twoCS2inv_pd,   ones_pd);
                    __m256d b_B              = _mm256_fmsub_pd( b_dv_dot,      b_A,            b_v_norm);
                    __m256d b_C              = _mm256_fmadd_pd( b_B,           twoCS2inv_pd,   ones_pd);
                    __m256d b_D              = _mm256_fmadd_pd( b_tmpB,        b_C,            b_tmpC);

                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index, D);
                    _mm256_storeu_pd(S->previous_particle_distributions+b_index, b_D);
            }
        }
    }
}



void collision_AVX4_u2(struct LBMarrays* S) {

    const int nBB = BLOCKSIZE_COL > S->nXYZ ? S->nXYZ : BLOCKSIZE_COL;
    const double tau_inv = 1/S->tau;
    const double cs2 = S->c_s*S->c_s;
    const double cs4 = cs2*cs2;


        __m256d ones_pd      = _mm256_set1_pd(1.0);
        __m256d CS2inv_pd    = _mm256_set1_pd(1.0/cs2);
        __m256d twoCS2inv_pd = _mm256_set1_pd(0.5/cs2);
        __m256d twoCS4inv_pd = _mm256_set1_pd(0.5/cs4);
        __m256d TAUinv_pd    = _mm256_set1_pd(tau_inv);
        __m256d omTAUinv_pd  = _mm256_set1_pd(1.0 - tau_inv);
        

        for(int nxyz = 0; nxyz < S->nXYZ; nxyz+=nBB) {
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
                    
                    
                    _mm256_storeu_pd(S->previous_particle_distributions+index,       D);
                    _mm256_storeu_pd(S->previous_particle_distributions+b_index,   b_D);
                    _mm256_storeu_pd(S->previous_particle_distributions+c_index,   c_D);
                    _mm256_storeu_pd(S->previous_particle_distributions+d_index,   d_D);
            }
        }
    }
}




// BLOCKING TEST ...
// choice between 256, 512, 1024 seems to be not relevant
