
#ifndef LBM_2d_H
#define LBM_2d_H


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

struct LBM{
    int Nx;
    int Ny;
    double* F;
    int* cylinder;
    const int* cxs;
    const int* cys;
    const double* weights;

    double *rho;
    double *ux;
    double *uy;
    double *Feq;
    double* temp_F;
    double* bndryF;
};


void perform_timestep(struct LBM* S);

double perform_measurement(int N, struct LBM* S, int time);

//struct LBM* lbm_initialize(int Nx, int Ny, const char* filename);

//void lbm_finalize(struct LBM* S);


#endif