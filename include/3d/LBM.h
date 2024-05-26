
#ifndef LBM_H_FILE
#define LBM_H_FILE

struct LBMarrays{
    int nX;
    int nY;
    int nZ;
    int nXY;
    int nXYZ;
    int direction_size;
    double* density_field;
    double* velocity_field;
    double* velocity_fieldX;
    double* velocity_fieldY;
    double* velocity_fieldZ;
    double* previous_particle_distributions;
    double* particle_distributions;
    int* reverse_indexes;
    const int* directions;
    const double* weights;
    double c_s;
    double c_s2;
    double c_s4;
    double tau;
    double tau_inv;
    double gamma_dot;
    int boundary_condition;
};


#endif