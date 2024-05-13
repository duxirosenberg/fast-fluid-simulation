
#include <stdio.h>
#include <stdlib.h>

#include "lbm2d_struct.h"
#include "utils2d.h"

#include "2D_drift.h"
#include "2D_set_reflective_bndry.h"
#include "2D_collision.h"
#include "2D_apply_bndry.h"



struct LBM* lbm_initialize(int Nx, int Ny, const char* filename) {
    struct LBM* S = (struct LBM*)malloc(sizeof(struct LBM));
    S->Nx = Nx;
    S->Ny = Ny;
    S->cxs = cxs;
    S->cys = cys;
    S->weights = weights;
    
    S->F = (double*)calloc(Nx * Ny * NL, sizeof(double));
    S->cylinder = (int*)calloc(Nx * Ny, sizeof(int));
    S->rho = (double*)calloc(Nx * Ny, sizeof(double));
    S->ux = (double*)calloc(Nx * Ny, sizeof(double));
    S->uy = (double*)calloc(Nx * Ny, sizeof(double));
    S->Feq = (double*)calloc(Nx * Ny * NL, sizeof(double));
    S->temp_F = (double*)malloc(Nx * Ny * NL * sizeof(double));
    S->bndryF = (double*)calloc(Nx * Ny * NL, sizeof(double));

    
    // Check if the filename is provided
    if (filename != NULL && strlen(filename) > 0) {
        // Read initial conditions from the provided file
        read_csv_to_array(filename, S->F, Nx * Ny * NL);
    } else {
        // Default initialization
        initialize_F(S->F, Nx, Ny);
    }
    initialize_cylinder(S->cylinder, Nx, Ny);

    return S;
}

void lbm_finalize(struct LBM* S) {
    free(S->F);
    free(S->cylinder);
    free(S->rho);
    free(S->ux);
    free(S->uy);
    free(S->Feq);
    free(S->temp_F);
    free(S->bndryF);
}

void perform_timestep(struct LBM* S) {
    drift(S);
    set_reflective_boundaries(S);
    collision(S);
    apply_boundary(S);
}

