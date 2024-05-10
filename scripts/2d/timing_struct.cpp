#include <iostream>
#include <vector>
#include <string>
#include <list>


#include "tsc_x86.h"
#include "utils2d.h"

#include "lbm2d_struct.h"

using namespace std;

#define CYCLES_REQUIRED 1e8
#define REP 5
#define EPS 1e-3



//Global Vars Drift
vector<f_structs> Funcs_drift;
vector<string> funcNames_drift;
vector<int> funcFlops_drift;
int numFuncs_drift = 0;
//Global Vars set reflective boundaries
vector<f_structs> Funcs_set_reflective_boundaries;
vector<string> funcNames_set_reflective_boundaries;
vector<int> funcFlops_set_reflective_boundaries;
int numFuncs_set_reflective_boundaries = 0;
//Global Vars collision
vector<f_structs> Funcs_collision;
vector<string> funcNames_collision;
int numFuncs_collision = 0;
vector<int> funcFlops_collision;
//Global Vars apply boundary
vector<f_structs> Funcs_apply_bndry;
vector<string> funcNames_apply_bndry;
vector<int> funcFlops_apply_bndry;
int numFuncs_apply_bndry = 0;



// Function to add new drift optimization
void add_function_drift(f_structs f, const char* name, int flops) {
    Funcs_drift.push_back(f);
    funcNames_drift.push_back(name);
    funcFlops_drift.push_back(flops);
    numFuncs_drift++;
}

// Function to add new reflective boundaries optimization
void add_function_set_reflective_boundaries(f_structs f, const char* name, int flops) {
    Funcs_set_reflective_boundaries.push_back(f);
    funcNames_set_reflective_boundaries.push_back(name);
    funcFlops_set_reflective_boundaries.push_back(flops);
    numFuncs_set_reflective_boundaries++;
}

// Function to add new collision optimization
void add_function_collision(f_structs f, const char* name, int flops) {
    Funcs_collision.push_back(f);
    funcNames_collision.push_back(name);
    funcFlops_collision.push_back(flops);
    numFuncs_collision++;
}

// Function to add new apply boundary optimization
void add_function_apply_boundary(f_structs f, const char* name, int flops) {
    Funcs_apply_bndry.push_back(f);
    funcNames_apply_bndry.push_back(name);
    funcFlops_apply_bndry.push_back(flops);
    numFuncs_apply_bndry++;
}



//Forward declaration of functions defined in lbm2d_baseline.c
extern "C" void register_drift_functions();
extern "C" void register_set_reflective_boundaries_functions();
extern "C" void register_collision_functions();
extern "C" void register_apply_boundary_functions();


extern "C" void lbm_baseline(double *F, int *cylinder, int Nx, int Ny, int n_steps);

typedef struct {
    double cycles;
    double performance;
} result;


result time_function(comp_func f, int flops, double* F, int* cylinder, int Nx, int Ny) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.  
    do{
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            //perform Nt steps of Simulation on initial conditions F_init, store the result in F
            f(F, cylinder, Nx, Ny,1);
        }
        end = stop_tsc(start);

        cycles = (double) end;
        multiplier = (CYCLES_REQUIRED) / cycles;
    } while (multiplier > 2);

    //reinitialize F
    initialize_F(F, Nx, Ny);
    // Perform REP runs and calculate the average number of cycles
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            //perform 1 step of Simulation
            f(F, cylinder, Nx, Ny, 1);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    double total_flops = flops * Nx * Ny;
    double performance  = total_flops / total_cycles;
    result res = {total_cycles, performance};
    return res; // Returns the average cycles over the repetitions
}

int main(int argc, char **argv) {

    register_drift_functions();
    register_set_reflective_boundaries_functions();
    register_collision_functions();
    register_apply_boundary_functions();


    cout << endl << "[1/] Starting program. " << numFuncs_drift + numFuncs_set_reflective_boundaries + numFuncs_collision + numFuncs_apply_bndry << " function(s) registered." << endl;
    cout << "Drift functions: " << numFuncs_drift << endl;
    cout << "Reflective boundaries functions: " << numFuncs_set_reflective_boundaries << endl;
    cout << "Collision functions: " << numFuncs_collision << endl;
    cout << "Apply boundary functions: " << numFuncs_apply_bndry << endl;

    /*-------------------------------------Benchmarking the Functions-------------------------------------------*/
    //iterate over various problem sizes, initialize F and cylinder, and run all functions. 
    //Nx = Ny
    int sizes[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    int numSizes = 10;
    //array of pointers to the various F and cylinder arrays
    struct LBM* lbmFuncs[numSizes];
    double* F_init[numSizes];
    double* F_baseline[numSizes];
    double* F_baseline2[numSizes];
    int* cylinder[numSizes];

    for (int i = 0; i < numSizes; i++) {
        F[i] = (double*)calloc(sizes[i]*sizes[i]*NL, sizeof(double));
        cylinder[i] = (int*)calloc(sizes[i]*sizes[i], sizeof(int));
        initialize_F(F[i], sizes[i], sizes[i]);
        initialize_cylinder(cylinder[i], sizes[i], sizes[i]);
    }


    //time function returns struct with total cycles and performance
    double results[numSizes][(numFuncs+1)*2];
    cout << "       Timing [1/" << numFuncs+1 << "]: Baseline" << endl;
    //baseline
    for (int i = 0; i < numSizes; i++) {
        int flops = 808;
        result res = time_function(lbm_baseline, flops, F[i], cylinder[i], sizes[i], sizes[i]);
        results[i][0] = res.cycles;
        results[i][1] = res.performance;
        cout << "       Timing [1/" << numFuncs+1 << "]: Baseline on average takes " << res.cycles << " cycles per timestep for a problem size of " << sizes[i] << "x" << sizes[i] << "." << endl;
        cout << "       Timing [1/" << numFuncs+1 << "]: Baseline has an average performance of " << res.performance << " FLOPS/cycle for a problem size of " << sizes[i] << "x" << sizes[i] << "." << endl;
    }
    //functions
    for (int j = 0; j < numFuncs; j++) {
        cout << "       Timing [" << j+1 << "/" << numFuncs+1 << "]: Function \"" << funcNames[j] << "\"" << endl;
        for (int i = 0; i < numSizes; i++) {
            result res = time_function(lbmFuncs[j], funcFlops[j], F[i], cylinder[i], sizes[i], sizes[i]);
            results[i][(j+1)*2] = res.cycles;
            results[i][(j+1)*2+1] = res.performance;
            cout << "       Timing [" << j+1 << "/" << numFuncs+1 << "]: Function \"" << funcNames[j] << "\" on average takes " << res.cycles << " cycles per timestep for a problem size of " << sizes[i] << "x" << sizes[i] << "." << endl;
            cout << "       Timing [" << j+1 << "/" << numFuncs+1 << "]: Function \"" << funcNames[j] << "\" has an average performance of " << res.performance << " FLOPS/cycle for a problem size of " << sizes[i] << "x" << sizes[i] << "." << endl;
        }
    }


    //ask whether user wants to store results in a csv 
    cout << "Do you want to store the results in a csv file? (y/n)" << endl;
    char answer;
    cin >> answer;
    if(answer == 'y' || answer == 'Y'){
        cout << "Enter the filename: " << endl;
        string filename;
        cin >> filename;
        cout << "Storing results in " << filename << endl;
        //store the results in a csv file
        FILE* file = fopen(filename.c_str(), "w");
        fprintf(file, "Problem size,");
        for (int j = 0; j < numFuncs+1; j++) {
            fprintf(file, "%s cycles,%s FLOPS/cycle,", j == 0 ? "Baseline" : funcNames[j-1].c_str(), j == 0 ? "Baseline" : funcNames[j-1].c_str());
        }
        fprintf(file, "\n");
        for (int i = 0; i < numSizes; i++) {
            fprintf(file, "%d,", sizes[i]);
            for (int j = 0; j < numFuncs+1; j++) {
                fprintf(file, "%f,%f,", results[i][j*2], results[i][j*2+1]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }

    //free memory
    for (int i = 0; i < numSizes; i++) {
        free(F[i]);
        free(cylinder[i]);
    }

    cout << "Program finished." << endl<<endl;
    
    return 0;
}
