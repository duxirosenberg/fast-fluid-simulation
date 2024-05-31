
#ifndef CMDLINE_LBM_COLLISION_H
#define CMDLINE_LBM_COLLISION_H
#include "timing.h"


// should be same as in collision.c
#define BLOCKSIZE_COL 256


// Flops: nX * nY * nZ * q * 29 + 2
// Intops: xN * nY * nZ * q * 23
// Doubles: 
            // prev_particle_dist: nXYZ*q;    write
            // velocity_field: nXYZ*3
            // denssity field nXYZ, 
            // weights q
//  Ints:   // dirs: 3xq
void collision_arrays(int nX, int nY, int nZ, int direction_size, double tau, double c_s,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights);

void collision_baseline(struct LBMarrays* S);


// basic strength reduction etc
void collision_1(struct LBMarrays* S);


// further strength reduction, 
// feq function inlined
// simplified triple loop of xyz to a simple loop
// avoid division
void collision_2(struct LBMarrays* S);


// reorder loops to the better memory dynamic variant
void collision_3(struct LBMarrays* S);


// cache optimized by blocking of loops 
// will only help for sufficiently large grids 
// (when the double Gridsize 2xN^3 doesnt fit into cache ... (usually at 32 or 64)
// Blocking only performs better starting fromm grid size 64 ..
//  currenlty chosen blocksize 256, exact blocksize has no "big" influece on performance
void collision_6(struct LBMarrays* S);


// implementation using SSA code design ...
void collision_SSA(struct LBMarrays* S); 
//unrolled once for(i+=2)(worse)
void collision_SSA_u(struct LBMarrays* S); 

// using the separated data structure for different velocity components
void collision_SSA2(struct LBMarrays* S); 
// unrolled once (worse)
void collision_SSA2_u(struct LBMarrays* S); 



// AVX2 VECTORISATION:
//  _u unrolled once  (+=8)
// _u2 unrolled twice (+=16)

// basic
void collision_AVX(struct LBMarrays* S);
void collision_AVX_u(struct LBMarrays* S);
void collision_AVX_u2(struct LBMarrays* S);

//avoid hadd() function in dot-product computation
void collision_AVX2(struct LBMarrays* S);

// use set_pd for loading and avoid complex dot product computation
void collision_AVX3(struct LBMarrays* S);
void collision_AVX3_u(struct LBMarrays* S);

// use new Datastructure with separated velocity components
void collision_AVX4(struct LBMarrays* S);
void collision_AVX4_u(struct LBMarrays* S);
void collision_AVX4_u2(struct LBMarrays* S);


// reordering seems to usually worsen, compiler seems to handle current
// instruction order pretty decently 
//  (even if it is definitely note the final instruction order 
//   that will be used)

///////////////////////////////////////////////////////////////////////////


static struct ops collision_baseline_flops(struct LBMarrays* S) {
    long val0 = S->nX * S->nY * S->nZ;
    long val = val0 * S->direction_size;
    struct ops ops = {
            2 + val * 29,
            val * 23,
            val0*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(3*sizeof(int) +sizeof(double)),
            val*(int)(sizeof(double))
    
    };
    return ops;
}


static struct ops collision_flops_3(struct LBMarrays* S) {
    long val = S->nXYZ * S->direction_size;
    struct ops ops = {
            5 + val * 22,
            val*5 + 3*S->direction_size,
            S->nXYZ*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(3*sizeof(int) +sizeof(double)),
            val*(int)(sizeof(double))
    
    };
    return ops;
}


static struct ops collision_flops_6(struct LBMarrays* S) {
    long val = S->nXYZ * S->direction_size;
    struct ops ops = {
            5 + val * 22,
            val*5 + (S->nXYZ/BLOCKSIZE_COL)*3*S->direction_size,
            S->nXYZ*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(3*sizeof(int) +sizeof(double)),
            val*(int)(sizeof(double))
    
    };
    return ops;
}
//same for SSA and AVX??




// string der ausgedruckt wird für rachel ...
// Collision_
// baseline, ,1, 2, 3, 4
// strength, memory, ssa, vecotrisation

static void register_collision_functions() {
    // add_collision_array_func(&collision_arrays, &collision_baseline_flops, "Collision - Arrays Bl");
    // add_collision_struct_func(&collision_baseline, &collision_baseline_flops, "Collision_0");  //baseline
    
 
    // for intel skylake, these 4 are optimal 4 steps ..
    // add_collision_struct_func(&collision_AVX4_u2, &collision_flops_6,    "Collision_4"); 
    // add_collision_struct_func(&collision_SSA2,    &collision_flops_6,    "Collision_3");
    // add_collision_struct_func(&collision_6,       &collision_flops_6,    "Collision_2");
    // add_collision_struct_func(&collision_3,       &collision_flops_3,    "Collision_1");

    add_collision_struct_func(&collision_3,      &collision_flops_3, "Collision 1");
    add_collision_struct_func(&collision_6,      &collision_flops_6, "Collision 2");

    // check
    add_collision_struct_func(&collision_SSA2,   &collision_flops_6, "Collision 3");
    //add_collision_struct_func(&collision_SSA2_u, &collision_flops_6, "Collision_SSA2_u");

    // check
    //add_collision_struct_func(&collision_AVX4,   &collision_flops_6, "Collision_AVX4");
    //add_collision_struct_func(&collision_AVX4_u, &collision_flops_6, "Collision_AVX4_u");
    add_collision_struct_func(&collision_AVX4_u2, &collision_flops_6, "Collision 4");
    

    //old data structure used unlikely to be better
    // add_collision_struct_func(&collision_SSA,    &collision_flops_6, "Collision_SSA");
    // add_collision_struct_func(&collision_SSA_u,  &collision_flops_6, "Collision_SSA_u");
    // add_collision_struct_func(&collision_AVX,    &collision_flops_6, "Collision_AVX");
    // add_collision_struct_func(&collision_AVX_u,  &collision_flops_6, "Collision_AVX_u");
    // add_collision_struct_func(&collision_AVX_u2, &collision_flops_6, "Collision_AVX_u2");   
    // add_collision_struct_func(&collision_AVX2,   &collision_flops_6, "Collision_AVX2");
    // add_collision_struct_func(&collision_AVX3,   &collision_flops_6, "Collision_AVX3");
    // add_collision_struct_func(&collision_AVX3_u, &collision_flops_6, "Collision_AVX3_u"); 


}

#endif //CMDLINE_LBM_COLLISION_H
