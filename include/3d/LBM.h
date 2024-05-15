
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
    double* previous_particle_distributions;
    double* particle_distributions;
    int* reverse_indexes;
    const int* directions;
    const double* weights;
    double c_s;
    double tau;
    double gamma_dot;
    int boundary_condition;
    //new fields for optimization of couette
    int* xDirections;
    int* yDirections;
    int* zDirections;
};


#endif