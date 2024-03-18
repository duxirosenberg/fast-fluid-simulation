#include <stdio.h>
#include <iostream>
#include <fstream>

#ifndef ASL_2024_OUTPUT_H
#define ASL_2024_OUTPUT_H

void output_velocity(int nX, int nY, double* velocity_field);

void output_lbm_data(std::string filename, bool header, int nX, int nY, int nZ, double* density_field, double* velocity_field);

void output_indices_file(int nX, int nY, int nZ);

#endif //ASL_2024_OUTPUT_H
