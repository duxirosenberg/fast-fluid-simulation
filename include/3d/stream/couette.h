
#ifndef CMDLINE_LBM_COUETTE_H
#define CMDLINE_LBM_COUETTE_H
#include "timing.h"

// Flops: 5 * (nX * nZ * q / 3)
// Intops: nX * nZ * q * (39 * nY + 32 / 3)
// Doubles: // particle_dist: nXYZ*q        write
            // prev_particle_dist: nXYZ*q;
            // weights q
//  Ints:   // dirs: 3xq
            // rev_indexes: q
void stream_couette_baseline(struct LBMarrays* S);
void stream_couette_code_motion(struct LBMarrays* S);
void stream_couette_loop_structure(struct LBMarrays* S);
void stream_couette_memcpy(struct LBMarrays* S);
void stream_couette_avx(struct LBMarrays* S);

void stream_couette_flattened(struct LBMarrays* S);
void stream_couette_flattened_avx(struct LBMarrays* S);



void stream_couette_arrays(int nX, int nY, int nZ, int direction_size, double c_s,
                           double* previous_particle_distributions,
                           double* particle_distributions,
                           const int* directions,
                           const double* weights,
                           const int* reverse_indexes);

void stream_couette_arrays_opt1(int nX, int nY, int nZ, int direction_size, double c_s,
                           double* previous_particle_distributions,
                           double* particle_distributions,
                           const int* directions,
                           const double* weights,
                           const int* reverse_indexes);




static struct ops stream_couette_baseline_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    struct ops ops = {
            (long)(5 * val / 3.),
            (long)(val * (39 * S->nY + 32 / 3.)),
            val*S->nY*2*(int)(sizeof(double)) + S->direction_size*(int)((3+1)*sizeof(int) +sizeof(double)),
            val*S->nY*(int)(sizeof(double))
    };
    return ops;
}


static struct ops stream_couette_codemotion_flops(struct LBMarrays* S) {
    long val = S->nX * S->nZ * S->direction_size;
    int q = S->direction_size;

    struct ops ops = {
        (long)(4 + q *(S->nZ * S->nX / 3.)), //flops
        (long)(q * (12*S->nXYZ+S->nXYZ-2*(S->nX*S->nZ*(1/3.)))), //iops
        (long)(2*q*(S->nXYZ)+27)*sizeof(double)+(6+3*q)*sizeof(int), //bytes read
        (long)(q*S->nXYZ)*sizeof(double) //bytes written
    };
    return ops;
}
 

static struct ops stream_couette_loop_structure_flops(struct LBMarrays* S) {
    int q = S->direction_size;
    struct ops ops = {
        (long)(4 + q *(2+ S->nZ * S->nX / 3.)), //flops
        (long)(q * (11+ S->nZ*(2+ S->nY*(1+S->nX*((7/3.)+2))))), //iops
        (long) (2*q*(S->nXYZ)+q)*sizeof(double)+(4*q)*sizeof(int), //bytes read
        (long)(q*S->nXYZ)*sizeof(double) //bytes written

    };
    return ops;
}


static struct ops stream_couette_memcpy_flops(struct LBMarrays* S) {
    int q = S->direction_size;
    int nZX = S->nZ * S->nX;
    int nXYZ = S->nX * S->nY * S->nZ;
    int nXY = S->nXY; int nZX = S->nZ*S->nX;
    int nX = S->nX;int nY = S->nY;int nZ = S->nZ;
    double t1;
    if(q==27){
        t1 = 3/27.;
    }else{//q==15
        t1 = 1/15.;
    }
    int t2 = q*((2/3.)*nZX + (3./q)*nXYZ + 2*t1*nZ*(nXY-nX) + (nZ/3.)*(nY-1/3)*nX);
    struct ops ops = {
        (long)(6 + q *(nZ * nX / 3)), //flops
        (long)(9 +q * ((1/3.)*nZ*6 + (1/3.)*(10+nX/2.) +4+4)), //iops
        (long) (1 +2*t2)*sizeof(double)+(10)*sizeof(int), //bytes read
        (long) t2*sizeof(double) //bytes written
    };
    return ops;
}







static void register_stream_couette_functions() {
    add_stream_couette_struct_func(&stream_couette_code_motion, &stream_couette_baseline_flops, "Couette 1");
    add_stream_couette_struct_func(&stream_couette_loop_structure, &stream_couette_loop_structure_flops, "Couette 2");
    add_stream_couette_struct_func(&stream_couette_memcpy, &stream_couette_memcpy_flops, "Couette 3");
    add_stream_couette_struct_func(&stream_couette_avx, &stream_couette_memcpy_flops, "Couette 4");
    
    //add_stream_couette_array_func(&stream_couette_arrays, &stream_couette_baseline_flops, "Stream Couette - Arrays Bl");
}

#endif //CMDLINE_LBM_COUETTE_H
