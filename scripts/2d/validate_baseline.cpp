#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>


#include "utils2d.h"


using namespace std;


extern "C" void lbm_baseline(double *F, int *cylinder, int Nx, int Ny, int n_steps);
// Stub version of add_function: For validate_baseline, this function does not do anything
extern "C" void add_function(comp_func f, const char* name, int flops){}


int main() {
    // Run the Python script to generate initial_conditions.csv and baseline.csv
    cout << endl << "Starting Validation of the C baseline:" << endl;
    cout << "[1/3] Running Python script (latticeboltzmann.py) to generate initial conditions and baseline"  << endl;
    system("python ../scripts/2d/latticeboltzmann.py");
    int Nx = 400, Ny=100; 
    
    // Allocate memory for the arrays
    double* F = (double*)calloc(Nx*Ny*NL, sizeof(double));
    double* result_baseline = (double*)calloc(Nx*Ny*NL, sizeof(double));
    int* cylinder = (int*)malloc(Nx*Ny * sizeof(int));
    
    //Read the initial conditions and the Python baseline
    read_csv_to_array("initial_conditions.csv", F, Nx * Ny * NL);
    read_csv_to_array("py_baseline.csv", result_baseline, Nx * Ny * NL);
    
    // Initialize the cylinder
    initialize_cylinder(cylinder, Nx, Ny);

    // Run the C baseline
    cout << "[2/3] Running the C baseline"<< endl;
    lbm_baseline(F, cylinder, Nx, Ny, Nt);

    // Compare the results to the Python Baseline
    cout << "[3/3] Comparing the results: " ;
    if (equal(result_baseline, F, Nx*Ny*NL)) {
        cout << "Results are equal" << endl;
    } else {
        cout << "Results are not equal" << endl;
    }

    //remove the files
    //remove("initial_conditions.csv");
    //remove("py_baseline.csv");

    cout << "Validation complete." << endl << endl;
    // Free the memory
    free(F);
    free(result_baseline);
    free(cylinder);
    return 0;
}
