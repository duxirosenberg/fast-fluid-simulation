#include <stdio.h>

#include <stdlib.h>
#include <math.h>

#define Nx 400    // resolution x-dir
#define Ny 100    // resolution y-dir
#define rho0 100    // average density
#define tau 0.6    // collision timescale
#define Nt 4000   // number of timesteps
#define plotRealTime 1 // switch on for plotting as the simulation goes along

// Lattice speeds / weights
#define NL 9
double cxs[] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
double cys[] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
double weights[] = {4.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36}; // sums to 1

// ---SAVING FUNCTIONS FOR CHECKING--------------//
void save_cylinder_to_csv(double *cylinder, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (fp != NULL) {
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                fprintf(fp, "%f,",  *(cylinder + i*Nx + j));
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    } else {
        printf("Error: Cannot open file for writing.\n");
    }
}

void save_F(double *F, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (fp != NULL) {
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                for (int k = 0; k < NL; k++) {
                    fprintf(fp, "%f,", *(F + i*Nx*NL + j*NL + k));
                }
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    } else {
        printf("Error: Cannot open file for writing.\n");
    }
}

// --------------------------------------------//



// Function to initialize F
void initialize_F(double *F) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < NL; k++) {
                *(F + i*Nx*NL + j*NL + k) = 1.0; // rho0 / NL;
            }
        }
    }

    // Add random noise
    srand(42);
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < NL; k++) {
                *(F + i*Nx*NL + j*NL + k) += 0.01 * ((double)rand() / RAND_MAX - 0.5);
            }
        }
    }

    // Add initial conditions for direction 3
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            *(F + i*Nx*NL + j*NL + 3) += 2 * (1 + 0.2 * cos(2 * M_PI * j / Nx * 4));
        }
    }

    // Normalize densities
    double rho;
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            rho = 0.0;
            for (int k = 0; k < NL; k++) {
                rho += *(F + i*Nx*NL + j*NL + k);
            }
            for (int k = 0; k < NL; k++) {
                *(F + i*Nx*NL + j*NL + k) *= rho0 / rho;
            }
        }
    }
}




int main() {
    // Allocate memory for F
    double *F = (double *)malloc(Nx * Ny * NL * sizeof(double));

    // Initialize F
    initialize_F(F);

    save_F(F, "initial_F1.csv");


    // Cylinder boundary
    double *cylinder = (double *)malloc(Nx * Ny * sizeof(double));
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            *(cylinder + i*Nx + j) = (pow(j - Nx/4, 2) + pow(i - Ny/2, 2) < pow(Ny/4, 2)) ? 1 : 0;
        }
    }

    save_cylinder_to_csv(cylinder, "cylinder2.csv");

    // Simulation Main Loop
    for (int it = 0; it < Nt; it++) {
        //printf("%d\n", it);

        // Drift
        for (int k = 0; k < NL; k++) {
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    int new_i = (i + (int)cys[k] + Ny) % Ny;
                    int new_j = (j + (int)cxs[k] + Nx) % Nx;
                    *(F + i*Nx*NL + j*NL + k) = *(F + new_i*Nx*NL + new_j*NL + k);
                }
            }
        }

        //save_F(F, "drift1.csv");

        // Set reflective boundaries
        double *bndryF = (double *)malloc(Nx * Ny * NL * sizeof(double));
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (*(cylinder + i*Nx + j) == 1) {
                    for (int k = 0; k < NL; k++) {
                        *(bndryF + i*Nx*NL + j*NL + k) = *(F + i*Nx*NL + j*NL + ((k+3)%NL));
                    }
                }
            }
        }

        //save_F(bndryF, "bndry1.csv");

        // Calculate fluid variables
        double *rho = (double *)malloc(Nx * Ny * sizeof(double));
        double *ux = (double *)malloc(Nx * Ny * sizeof(double));
        double *uy = (double *)malloc(Nx * Ny * sizeof(double));
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                double temp_rho = 0.0;
                double temp_ux = 0.0;
                double temp_uy = 0.0;
                for (int k = 0; k < NL; k++) {
                    temp_rho += *(F + i*Nx*NL + j*NL + k);
                    temp_ux += *(F + i*Nx*NL + j*NL + k) * cxs[k];
                    temp_uy += *(F + i*Nx*NL + j*NL + k) * cys[k];
                }
                *(rho + i*Nx + j) = temp_rho;
                *(ux + i*Nx + j) = temp_ux / temp_rho;
                *(uy + i*Nx + j) = temp_uy / temp_rho;
            }
        }

        // Apply Collision
        double *Feq = (double *)malloc(Nx * Ny * NL * sizeof(double));
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                for (int k = 0; k < NL; k++) {
                    double cu = cxs[k] * (*(ux + i*Nx + j)) + cys[k] * (*(uy + i*Nx + j));
                    double uu = (*(ux + i*Nx + j)) * (*(ux + i*Nx + j)) + (*(uy + i*Nx + j)) * (*(uy + i*Nx + j));
                    *(Feq + i*Nx*NL + j*NL + k) = *(rho + i*Nx + j) * weights[k] * (1 + 3 * cu + 9 * cu * cu / 2 - 3 * uu / 2);
                }
            }
        }

        //save_F(Feq, "feq1.csv");

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                for (int k = 0; k < NL; k++) {
                    *(F + i*Nx*NL + j*NL + k) -= (1.0 / tau) * (*(F + i*Nx*NL + j*NL + k) - *(Feq + i*Nx*NL + j*NL + k));
                    //printf("%f\n", *(F + i*Nx*NL + j*NL + k));
                }
            }
        }
        

        // Apply boundary
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (*(cylinder + i*Nx + j) == 1) {
                    for (int k = 0; k < NL; k++) {
                        *(F + i*Nx*NL + j*NL + k) = *(bndryF + i*Nx*NL + j*NL + k);
                        //printf("%f\n", *(F + i*Nx*NL + j*NL + k));

                    }
                }
            }
        }
        //save_F(F, "fend1.csv");


        

        

        // Free memory
        free(bndryF);
         free(rho);
        free(ux);
        free(uy);
        free(Feq);
    }


    


    // Free memory
    free(F);
    free(cylinder);


    return 0;

}