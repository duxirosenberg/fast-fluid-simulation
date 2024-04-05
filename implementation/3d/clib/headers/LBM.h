
#ifndef LBM_H_FILE
#define LBM_H_FILE

void initialize(int nX, int nY, int nZ, int direction_size,
                double* density_field,
                double* velocity_field,
                double* previous_particle_distributions,
                double* particle_distributions,
                const int* directions,
                const double* weights,
                int* reverse_indexes
);

void perform_timestep(int nX, int nY, int nZ, int direction_size, int time, double tau, double gamma_dot, double c_s, int boundary_condition,
                      double* density_field,
                      double* velocity_field,
                      double* previous_particle_distributions,
                      double* particle_distributions,
                      const int* directions,
                      const double* weights,
                      int* reverse_indexes);

#endif