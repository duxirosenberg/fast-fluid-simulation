
#ifndef CMDLINE_LBM_COLLISION_H
#define CMDLINE_LBM_COLLISION_H
#include "timing.h"

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

// reduce ops by c_s powers and index precomputation

// basic strengthening via multiplication reduction etc ...
void collision_1(struct LBMarrays* S);

// simplified triple loop of xyz to a simple loop
void collision_2(struct LBMarrays* S);


// set as current baseline ...
//  optimized by simply switching the two loops
void collision_3(struct LBMarrays* S);

// only performs better starting fromm grid size 64 ..
// attempts at blocking all performing worse in most cases...
void collision_4(struct LBMarrays* S);
void collision_5(struct LBMarrays* S);
void collision_6(struct LBMarrays* S);
void collision_7(struct LBMarrays* S);
void collision_8(struct LBMarrays* S);
void collision_9(struct LBMarrays* S);


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


static struct ops collision_flops_2(struct LBMarrays* S) {
    long val0 = S->nX * S->nY * S->nZ;
    long val = val0 * S->direction_size;
    struct ops ops = {
            2 + val * 23,
            val * 8 + 5*S->direction_size,
            val0*(S->direction_size+3+1)*(int)(sizeof(double)) + S->direction_size*(int)(3*sizeof(int) +sizeof(double)),
            val*(int)(sizeof(double))
    
    };
    return ops;
}


static void register_collision_functions() {
    // add_collision_array_func(&collision_arrays, &collision_baseline_flops, "Collision - Arrays Bl");
    // add_collision_struct_func(&collision_baseline, &collision_baseline_flops, "Collision_0");  //baseline
    // add_collision_struct_func(&collision_2, &collision_flops_2, "Collision_2");   // basic strength increase, function inline etc
    // add_collision_struct_func(&collision_3, &collision_flops_2, "Collision_3");   // loop reordering for more efficient memory dynamic   
    // add_collision_struct_func(&collision_6, &collision_flops_2, "Collision_6");   // blocking of loops, will only help for sufficiently large gridss 1024



    // add_collision_struct_func(&collision_1, &collision_baseline_flops, "Collision_1");
    // add_collision_struct_func(&collision_4, &collision_baseline_flops, "Collision_4"); Blocking with 256
    // add_collision_struct_func(&collision_5, &collision_baseline_flops, "Collision_5"); Blocking with 512
    // add_collision_struct_func(&collision_5, &collision_baseline_flops, "Collision_7"); Blocking with 2048

    add_collision_struct_func(&collision_8, &collision_baseline_flops, "Collision_8"); // vectorisation works but non efficient
    add_collision_struct_func(&collision_9, &collision_baseline_flops, "Collision_9");  // doesnt work yet sadly ...
}

#endif //CMDLINE_LBM_COLLISION_H
