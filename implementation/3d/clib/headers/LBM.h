
#ifndef LBM_H_FILE
#define LBM_H_FILE


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#ifdef LBM_STRUCT


struct LBM{
    int nX;
    int nY;
    int nZ;
    int direction_size;
    double* density_field;
    double* velocity_field;
    double* previous_particle_distributions;
    double* particle_distributions;
    int* reverse_indexes;
    const int* directions;
    const double* weights;
    double c_s;
    double tau;
    double gamma_dot;
    const int boundary_condition;
};


void perform_timestep(struct LBM* S, int time);

double perform_Measurement(int N, struct LBM* S, int time);

#else



void perform_timestep(int nX, int nY, int nZ, int direction_size, int time, double tau, double gamma_dot, double c_s, int boundary_condition,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights,
                      int* reverse_indexes);


double perform_Measurement(int N, int nX, int nY, int nZ, int direction_size, 
                      int  time, double tau, double gamma_dot, double c_s, int boundary_condition,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights,
                      int* reverse_indexes);

#endif

void initialize(int nX, int nY, int nZ, int direction_size,
                double* density_field,
                double* velocity_field,
                double* previous_particle_distributions,
                double* particle_distributions,
                const int* directions,
                const double* weights,
                int* reverse_indexes);

#endif