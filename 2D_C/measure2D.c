
//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
//#include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


#ifdef __x86_64__
#include "tsc_x86.h"
#endif

#ifdef TEST
#include <assert.h>
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2e9
#define CALIBRATE


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




int compute() {
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


/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(int n) {
    int i, num_runs;
    myInt64 cycles = 0;
    myInt64 start;
    num_runs = NUM_RUNS;

    /*
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        for (i = 0; i < num_runs; ++i) {
            start = start_tsc();
            compute(n);
            cycles += stop_tsc(start);
        }

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    cycles = 0;
    for (i = 0; i < num_runs; ++i) {
        start = start_tsc();
        compute();
        cycles += stop_tsc(start);
    }

    cycles = cycles/num_runs;
    return (double) cycles;
}
#endif

double c_clock(int n) {
    int i, num_runs;
    double cycles = 0;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        for (i = 0; i < num_runs; ++i) {
            start = clock();
            compute();
            end = clock();
            cycles += (double)(end - start);
        }

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of CLOCKS_PER_SEC
        if(cycles >= CYCLES_REQUIRED/(FREQUENCY/CLOCKS_PER_SEC)) break;

        num_runs *= 2;
    }
#endif
    cycles = 0;
    for(i=0; i<num_runs; ++i) {
        start = clock();
        compute();
        end = clock();
        cycles += (double)(end - start);
    }
    return (cycles)/num_runs;
}

#ifndef WIN32
double timeofday(int n) {
    int i, num_runs;
    double cycles = 0;
    struct timeval start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        for (i = 0; i < num_runs; ++i) {
            gettimeofday(&start, NULL);
            compute();
            gettimeofday(&end, NULL);
            cycles += (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)*FREQUENCY;
        }

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    double seconds = 0;
    for(i=0; i < num_runs; ++i) {
        gettimeofday(&start, NULL);
        compute();
        gettimeofday(&end, NULL);
        seconds += (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6);
    }

    return seconds/num_runs;
}

#else

double gettickcount(int n) {
    int i, num_runs;
    double cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {

        for (i = 0; i < num_runs; ++i) {
            start = (double)GetTickCount();
            compute();
            end = (double)GetTickCount();
            cycles += (double)(end - start);
        }

        cycles = cycles*FREQUENCY/1e3; // end-start provides a measurement in the order of milliseconds

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    cycles = 0;
    for(i=0; i < num_runs; ++i) {

            start = (double)GetTickCount();
            compute();
            end = (double)GetTickCount();
            cycles += (double)(end - start);
    }

    return cycles/num_runs;
}

double queryperfcounter(int n, LARGE_INTEGER f) {
    int i, num_runs;
    double cycles = 0;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        for (i = 0; i < num_runs; ++i) {
            QueryPerformanceCounter(&start);
            compute();
            QueryPerformanceCounter(&end);
            cycles += (double)(end.QuadPart - start.QuadPart);
        }

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of f
        if(cycles >= CYCLES_REQUIRED/(FREQUENCY/f.QuadPart)) break;

        num_runs *= 2;
    }
#endif

    cycles = 0;
    for(i=0; i < num_runs; ++i) {
        QueryPerformanceCounter(&start);
        compute();
        QueryPerformanceCounter(&end);
        cycles += (double)(end.QuadPart - start.QuadPart);
    }

    return cycles/num_runs;
}

#endif

int main(int argc, char **argv) {
    if (argc!=2) {printf("usage: FW <n>\n"); return -1;}
    int n = atoi(argv[1]);
    printf("n=%d \n",n);

    //double* A = (double *)malloc(n*n*sizeof(double));
    //double* L = (double *)malloc(n*n*sizeof(double));
    //double* U = (double *)malloc(n*n*sizeof(double));

    //fill_matrix(A, n);

#ifdef __x86_64__
    double r = rdtsc(n);
    printf("RDTSC instruction:\n %lf cycles measured => %lf seconds, assuming frequency is %lf MHz. (change in source file if different)\n\n", r, r/(FREQUENCY), (FREQUENCY)/1e6);
#endif

    double c = c_clock(n);
    printf("C clock() function:\n %lf cycles measured. On some systems, this number seems to be actually computed "
           "from a timer in seconds then transformed into clock ticks using the variable CLOCKS_PER_SEC. Unfortunately, "
           "it appears that CLOCKS_PER_SEC is sometimes set improperly. (According to this variable, your computer should "
           "be running at %lf MHz). In any case, dividing by this value should give a correct timing: %lf seconds. \n\n",c, (double) CLOCKS_PER_SEC/1e6, c/CLOCKS_PER_SEC);


#ifndef WIN32
    double t = timeofday(n);
    printf("C gettimeofday() function:\n %lf seconds measured\n\n",t);
#else
    LARGE_INTEGER f;
    double t = gettickcount(n);
    printf("Windows getTickCount() function:\n %lf milliseconds measured\n\n",t);

    QueryPerformanceFrequency((LARGE_INTEGER *)&f);

    double p = queryperfcounter(n, f);
    printf("Windows QueryPerformanceCounter() function:\n %lf cycles measured => %lf seconds, with reported CPU frequency %lf MHz\n\n",p,p/f.QuadPart,(double)f.QuadPart/1000);
#endif

    return 0;
}









