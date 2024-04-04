#include <iomanip>
#include "output.hpp"


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
//Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//RapidJSON files.
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"
#include "csv.h"


#include <chrono>
#include <thread>


inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


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

void output_indices_file(std::string filename,int nX, int nY, int nZ) {
    std::ofstream output(filename);
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


void check_solution(std::string soluFolder,std::string testFolder, int time, int nX, int nY, int nZ, double* density_field, double* velocity_field){

        const std::string  soluFile = soluFolder+std::to_string(time)+".csv";
        const std::string  testFile = testFolder+std::to_string(time)+".csv";
        std::cout << "Comparing:  " << soluFile << "  to  " << testFile << "  at time  " << time << std::endl;

        bool equality = true; 
        if(file_exists(soluFile) /*&& file_exists(testFile)*/ ) {
                
        		io::CSVReader<4> inSolu(soluFile);
        		io::CSVReader<4> inTest(testFile);
        		inTest.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
        		inSolu.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
        		double density,u_x,u_y,u_z;
        		double densityS,u_xS,u_yS,u_zS;
        		for(int x = 0; x < nX; x++) {
        			for(int y = 0; y < nY; y++) {
        				for(int z = 0; z < nZ; z++) {
                            int index = (z * nX * nY) + (y * nX) + x;


                            inTest.read_row(density,u_x,u_y,u_z);
                            inSolu.read_row(densityS,u_xS,u_yS,u_zS);

                            equality =  equality
                                                && densityS == density
                                                && u_xS == u_x
                                                && u_yS == u_y
                                                && u_zS == u_z;
        				}
        			}
        		}
        	} else {
        		std::cout << "The solution File: " << soluFile << "does not exist" << '\n';
                return;
            }

            if(equality){
                std::cout<< "The Testresult at time " << time << " coincides with the Solution" << std::endl;
            }else{
                std::cout<< "The Testresult at time " << time << " does not coincide with the Solution" << std::endl;
            }

            return;
}



//doesnt work like this
// void check_solution(std::string soluFolder,, int time, int nX, int nY, int nZ, double* density_field, double* velocity_field){

//         const std::string  soluFile = soluFolder+std::to_string(time)+".csv";

//         bool equality = true; 
//         if(file_exists(soluFile)) {
//         std::cout << soluFile << std::endl;
//         		std::cout << "Checking Solution for Time: " << time << std::endl; // '\n';
//         		io::CSVReader<4> in(soluFile);
//         		in.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
//         		double density,u_x,u_y,u_z;
//         		for(int x = 0; x < nX; x++) {
//         			for(int y = 0; y < nY; y++) {
//         				for(int z = 0; z < nZ; z++) {
//                             int index = (z * nX * nY) + (y * nX) + x;
//                             in.read_row(density,u_x,u_y,u_z);

//                             equality =  equality
//                                                 && (density_field[index] == density)
//                                                 && (velocity_field[3 * index] == u_x)
//                                                 && (velocity_field[3 * index + 1] == u_y)
//                                                 && (velocity_field[3 * index + 2] == u_z);
//                             std::cout 
//                                                 << " " <<  density_field[index]          << " " << density
//                                                 << " " <<  velocity_field[3 * index]     << " " << u_x
//                                                 << " " <<  velocity_field[3 * index + 1] << " " << u_y
//                                                 << " " <<  velocity_field[3 * index + 2] << " " << u_z << std::endl;
//         				}
//         			}
//         		}
//         	} else {
//         		std::cout << "The solution File: " << soluFile << "does not exist" << '\n';
//                 return;
//             }

//             if(equality){
//                 std::cout<< "The Testresult at time " << time << "coincides with the Solution" << std::endl;
//             }else{
//                 std::cout<< "The Testresult at time " << time << "does not coincide with the Solution" << std::endl;
//             }

//             return;
// }