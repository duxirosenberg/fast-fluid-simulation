#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
//Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "output.hpp"
//RapidJSON files.
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"
#include "csv.h"

extern "C" {
    #include "LBM.h"
}

using namespace rapidjson;

inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

const int D2Q9_DIRECTIONS[] = { 1,0,0, 0,1,0, -1,0,0, 0,-1,0, 1,1,0, -1,1,0, -1,-1,0, 1,-1,0, 0,0,0};
const double D2Q9_WEIGHTS[] = {1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,4.0/9.0};

const int D3Q15_DIRECTIONS[] = { 0,0,0, 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1, 1,1,1, -1,-1,-1, 1,1,-1, -1,-1,1, 1,-1,1, -1,1,-1, -1,1,1, 1,-1,-1 };
const double D3Q15_WEIGHTS[] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0  };

const int D3Q27_DIRECTIONS[] = { 0,0,0,  1,0,0, -1,0,0,0,1,0, 0,-1,0, 0,0,1,
                                 0,0,-1, 1,1,0, -1,-1,0,1,0,1,-1,0,-1, 0,1,1,
                                 0,-1,-1, 1,-1,0, -1,1,0,1,0,-1, -1,0,1,0,1,-1,
                                 0,-1,1, 1,1,1, -1,-1,-1,1,1,-1, -1,-1,1, 1,-1,1,
                                 -1,1,-1, -1,1,1, 1,-1,-1};
const double D3Q27_WEIGHTS[] = {8.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0, 1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0,1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};


int main(int argc, char** argv) {
  if(!file_exists("options.json")) {
    std::cout << "Please ensure that options.json exists. If not, it can be obtained from root directory of GitHub repo." << '\n';
    return -1;
  }

	std::ifstream t("options.json");
	std::string str((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());
	rapidjson::Document d;
	d.Parse(str.c_str());
	

	std::string testFolder = d["outputFolder"].GetString();
	std::string soluFolder = d["solutionFolder"].GetString();

	#ifdef LBM_STRUCT
	testFolder.pop_back();
	testFolder = testFolder + "Struct/";
	#endif


	std::cout << "Do you want to clean the previous run? (1 - Yes, 0 - No): ";
	int choice=1;
	// std::cin >> choice;
	if(choice == 1) {
		std::string remove_command = "rm -rf ./" + testFolder;
		std::string mkdir_command = "mkdir ./" + testFolder;
		std::string copy_command = "cp options.json ./" + testFolder;

		system( remove_command.c_str());
		system( mkdir_command.c_str());
		system( copy_command.c_str());
	}

	Value& save_every_value = d["save_every"];
	int save_every = save_every_value.GetInt();
	std::cout << "Save every: " << save_every << '\n';
	auto grid_size = d["grid_size"].GetArray();
	int nX = grid_size[0].GetInt();
	int nY = grid_size[1].GetInt();
	int nZ = grid_size[2].GetInt();
	std::cout << "Grid size: " << nX << "x" << nY << "x" << nZ << '\n';
	Value& m_c_s = d["c_s"];
	double c_s = m_c_s.GetDouble();
	std::cout << "c_s (Speed of sound): " << c_s << '\n';
	Value& m_tau = d["tau"];
	double tau = m_tau.GetDouble();
	std::cout << "Tau: " << tau << '\n';
	if(tau < 0.5) {
	    //Section 3.5.5.1 with delta t=1.
	    std::cout << "Error: Tau must be greater than 0.5 for numerical stability." << '\n';
	    return -1;
	}
	Value& m_velocity_set = d["velocity_set"];
	std::string velocity_set = m_velocity_set.GetString();
	Value& m_boundary_conditions = d["boundary_conditions"];
	std::string boundary_conditions = m_boundary_conditions.GetString();
    int direction_size;

    if(velocity_set == "D3Q15") {
        direction_size = 15;
    } else if(velocity_set == "D3Q27") {
        direction_size = 27;
    } else if(velocity_set == "D2Q9") {
        direction_size = 9;
    } else {
        std::cout << "Error: Please specify a valid velocity set such as D3Q15,D3Q27 or D2Q9." << '\n';
        return -1;

    }

	std::cout << "Velocity set: " << velocity_set << '\n';

    int boundary_condition;
    if(boundary_conditions == "periodic") {
        boundary_condition = 1;
    } else if (boundary_conditions == "couette"){
        boundary_condition = 2;
    } else if(boundary_conditions == "lees_edwards") {
        boundary_condition = 3;
    } else {
        std::cout << "Errors: boundary_conditions in options.json can either be: periodic, Couette (D2Q9 only) or lees_edwards (Lees-Edwards Shear, Please see research paper by Alexander Wagner)";
        return -1;
    }
	std::cout << "Boundary conditions: " << boundary_conditions << '\n';
	    Value& m_gamma_dot = d["gamma_dot"];
	    double gamma_dot = m_gamma_dot.GetDouble();
	std::cout << "Shear rate (gamma_dot): " << gamma_dot << '\n';
	if(nZ != 1 && velocity_set == "D2Q9") {
	    std::cout << "Warning: NZ=1 for D2Q9.";
	    return -1;
	}
    Value& m_n_steps = d["n_steps"];
    int n_steps = m_n_steps.GetInt();

    std::cout << "Initializing Arrays" << std::endl;



    int box_flatten_length = nX * nY * nZ;
    int distributions_flatten_length = box_flatten_length * direction_size;

    double* density_field = (double*) malloc(box_flatten_length * sizeof(double));
    double* velocity_field = (double*) malloc(3 * box_flatten_length * sizeof(double));
    double* previous_particle_distributions = (double*) malloc(distributions_flatten_length * sizeof(double));
    double* particle_distributions = (double*) malloc(distributions_flatten_length * sizeof(double));
    int* reverse_indexes = (int*) malloc(direction_size * sizeof(int));
    const int* directions;
    const double* weights;

    if(direction_size == 15) {
        directions = &D3Q15_DIRECTIONS[0];
        weights = &D3Q15_WEIGHTS[0];
    } else if(direction_size == 27) {
        directions = &D3Q27_DIRECTIONS[0];
        weights = &D3Q27_WEIGHTS[0];
    } else if(direction_size == 9) {
        // 2D case
        directions = &D2Q9_DIRECTIONS[0];
        weights = &D2Q9_WEIGHTS[0];
    }
    initialize(nX, nY, nZ, direction_size, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);

    std::cout << "Finished Initializing Arrays" << std::endl;
	for(int i = 0; i < argc; i++) {
		if(std::string(argv[i]) == "generate_ic") {
			output_lbm_data("ic.csv", true, nX, nY, nZ, density_field, velocity_field);
			std::cout << "Generated ic.csv" << '\n';
			return 0;
		}
	}


	if(file_exists("ic.csv")) {
		std::cout << "Loading initial conditions" << '\n';
		io::CSVReader<4> in("ic.csv");
		in.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
		double density,u_x,u_y,u_z;
		for(int x = 0; x < nX; x++) {
			for(int y = 0; y < nY; y++) {
				for(int z = 0; z < nZ; z++) {
                    int index = (z * nX * nY) + (y * nX) + x;
					in.read_row(density,u_x,u_y,u_z);
                    density_field[index] = density;
                    velocity_field[3 * index] = u_x;
                    velocity_field[3 * index + 1] = u_y;
                    velocity_field[3 * index + 2] = u_z;
				}
			}
		}
		std::cout << "Loaded initial conditions" << '\n';
	} else {
		std::cout << "Using default of p=1 for all x,y,z and u(x,t=0)=0 for all x,y,z. (Steady state)" << '\n';
		std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<NX,j<NY and k<NZ" << '\n';
	}
	//Equation 3.5 with delta t = 1, LBM Principles and Practice book.
	double viscosity = c_s * c_s * (tau - 0.5);
	std::cout << "Kinematic shear viscosity: " << viscosity << '\n';
    //Equation 4.4.49
    std::cout << "For velocity set D2Q9,D3Q15 and D3Q27, |u_max|<0.577\n";


	int scale = 1;
	int runs = n_steps * scale * scale * scale;


	#ifdef LBM_STRUCT
	LBM solver{	nX, nY, nZ, direction_size, density_field, velocity_field, 
				previous_particle_distributions, particle_distributions, reverse_indexes, directions, weights, 
				c_s, tau, gamma_dot, boundary_condition};

	LBM* S =  &solver;
	std::cout << "solver allocated" << std::endl;
	#endif


	output_lbm_data(testFolder + "0.csv", true, nX, nY, nZ, density_field, velocity_field);
	output_indices_file(testFolder + "indices.csv", nX, nY, nZ);

	double measure;
	int time = 1; // time
	int* a = &time;


	__sync_synchronize(); 

		measure =
		#ifdef LBM_STRUCT 
		perform_Measurement(runs, S, time );
		#else
		perform_Measurement(runs, nX, nY, nZ, direction_size, time, tau, gamma_dot, c_s, boundary_condition, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);
		#endif
		time += runs-1;

	 __sync_synchronize(); 
	std::cout << std::endl;

	output_lbm_data(testFolder + std::to_string(time) + ".csv", true, nX, nY, nZ, density_field, velocity_field);	
	check_solution(soluFolder,testFolder, time, nX,nY,nZ, density_field, velocity_field);
    std::cout << "The simulation was performed with an average speed of  "  << measure << " mili-seconds per iteration" << std::endl;;

	return 0;
}







    // for(int i = 1; i <= runs; i++) {		
	// 	#ifdef LBM_STRUCT
	// 	perform_timestep(S, i);
	// 	#else
	// 	perform_timestep(nX, nY, nZ, direction_size, i, tau, gamma_dot, c_s, boundary_condition, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);
	// 	#endif

	// 	if((i) % save_every == 0) {
    //         double percentage = (double) (i) / (double) (runs) * 100.0;
    //         std::cout << "Saving data - " << (i) << "/" << runs << " (" << percentage << "%)" << '\n';
	// 		output_lbm_data(testFolder + std::to_string(i) + ".csv", true, nX, nY, nZ, density_field, velocity_field);
			
    //         //output_velocity(nX, nY, velocity_field);
    //     }
	// }
	// check_solution(soluFolder,testFolder, runs, nX,nY,nZ, density_field, velocity_field);