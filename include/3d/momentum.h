
#ifndef CMDLINE_LBM_MOMENTUM_H
#define CMDLINE_LBM_MOMENTUM_H
#include "timing.h"

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
// Doubles: // particle_dist: nXYZ*q;
            // velocity_field: nXYZ*3  write
            // denssity field nXYZ,    write  
// Ints    // dirs: 3xq
void momentum_baseline(struct LBMarrays* S);


void momentumO1(struct LBMarrays* S);
void momentumO2(struct LBMarrays* S);
void momentumO25(struct LBMarrays* S);


//void momentumO21(struct LBMarrays* S);
//void momentumO22(struct LBMarrays* S);
//void momentumO3(struct LBMarrays* S);
//void momentum_arraysO2(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions);
//void momentumO11(struct LBMarrays* S);
//void momentumO12(struct LBMarrays* S);
//void momentumO23(struct LBMarrays* S);
//void momentumO24(struct LBMarrays* S);
//void momentumO26(struct LBMarrays* S);

void momentum_arrays(int nX, int nY, int nZ, int direction_size, double* density_field, double* velocity_field, double* particle_distributions, const int* directions);

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: nX * nY * nZ * (10 + 14 * q)
static struct ops momentum_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (3 + 4 * S->direction_size),
            val * (10 + 14 * S->direction_size),
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}

// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: 3 + nX * nY (1 + nZ * (5 + q))
static struct ops momentum_O11_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (3 + 4 * S->direction_size),
            3 + val * (5 + S->direction_size) + S->nX * S->nY,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}
// Flops: nX * nY * nZ * (3 + 4 * q)
// Intops: 3 + nX * nY (1 + nZ * (3 + 4 * q))
static struct ops momentum_O12_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (3 + 4 * S->direction_size),
            3 + val * (3 + 4*S->direction_size) + S->nX * S->nY,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}

// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 3 + nX * nY (1 + nZ * (3 + 4 * q))
static struct ops momentum_O1_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (4 + 4 * S->direction_size),
            3 + val * (3 + 4*S->direction_size) + S->nX * S->nY,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}

// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 2 + q + 4 * q * nX * nY * nZ
static struct ops momentum_O2_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (4 + 4 * S->direction_size),
            2 + S->direction_size + S->direction_size * val * 4,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}




// Flops: nX * nY * nZ * (4 + 4 * q)
// Intops: 3 + nZ(1 + nY(1 + nX(4 + 5q))) 
static struct ops momentum_O21_flops(struct LBMarrays* S) {
    long val = S->nX * S->nY * S->nZ;
    struct ops ops = {
            val * (4 + 4 * S->direction_size),
            3 + val * (4 + 5*S->direction_size) + S->nZ * S->nY + S->nZ,
            val*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(sizeof(int)),
            val*4*(int)(sizeof(double))
    };
    return ops;
}

static void register_momentum_functions() {
    add_momentum_array_func(&momentum_arrays, &momentum_baseline_flops, "Momentum - Arrays Bl");
    add_momentum_struct_func(&momentum_baseline, &momentum_baseline_flops, "Momentum - Structs Bl");
    
    add_momentum_struct_func(&momentumO1, &momentum_O1_flops, "Momentum - 01 - precalculations and strength reduction"); //best in optimization step 1
    add_momentum_struct_func(&momentumO2, &momentum_O2_flops, "Momentum - 02 - best loop reorder etc"); 
    add_momentum_struct_func(&momentumO25, &momentum_O2_flops, "Momentum - 025 - best loop order w unrolling");


    //add_momentum_array_func(&momentum_arraysO2, &momentum_baseline_flops, "Momentum - 02 arrays");
    //add_momentum_struct_func(&momentumO3, &momentum_baseline_flops, "Momentum - 03 - AVX");
    //add_momentum_struct_func(&momentumO11, &momentum_O11_flops, "Momentum - O11");
    //add_momentum_struct_func(&momentumO12, &momentum_O12_flops, "Momentum - O12 ");
    //add_momentum_struct_func(&momentumO23, &momentum_O21_flops, "Momentum - 023 - more reordering");
    //add_momentum_struct_func(&momentumO21, &momentum_O21_flops, "Momentum - 021 - unroll x4");
    //add_momentum_struct_func(&momentumO22, &momentum_O2_flops, "Momentum - 022 - unroll x8");
    //add_momentum_struct_func(&momentumO24, &momentum_O2_flops, "Momentum - 04 - new loop order wo unroll");


}

#endif //CMDLINE_LBM_MOMENTUM_H
