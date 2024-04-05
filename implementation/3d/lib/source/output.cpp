#include <iomanip>
#include "output.hpp"


// TODO johannes: Seems like this method was only designed for 2D velocities?
void output_velocity(int nX, int nY, double* velocity_field) {
    int z_index = 0;
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            int index = (y * nX) + x;
            std::cout << velocity_field[3 * index] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void output_indices_file(int nX, int nY, int nZ) {
    std::ofstream output("output/indices.csv");
    output << "x,y,z" << '\n';
    for(int i = 0; i < nX; i++) {
        for(int j = 0; j < nY; j++) {
            for(int k = 0; k < nZ; k++) {
                output << i << "," << j << "," << k << '\n';
            }
        }
    }
    std::cout << std::endl;
    output.close();
}

void output_lbm_data(std::string filename, bool header, int nX, int nY, int nZ, double* density_field, double* velocity_field) {
    std::ofstream output_stream;
    output_stream.open (filename, std::ofstream::out | std::ofstream::app);
    if(header) {
        output_stream << "p,u_x,u_y,u_z" << '\n';
    }
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            for(int z = 0; z < nZ; z++) {
                int index = (z * nX * nY) + (y * nX) + x;
                output_stream << density_field[index] << "," <<
                              velocity_field[3 * index] << "," <<
                              velocity_field[3 * index + 1] << "," <<
                              velocity_field[3 * index + 2] << '\n';
            }
        }
    }
    output_stream.close();
}