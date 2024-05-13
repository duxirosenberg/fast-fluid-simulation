#include <iostream>
#include <vector>
#include <string>
#include <list>
#include "tsc_x86.h"
#include "utils2d.h"

using namespace std;

#define CYCLES_REQUIRED 1e8
#define REP 5
#define EPS 1e-3

// Global vars, used to keep track of LBM simulation functions
vector<comp_func> lbmFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

// Function to add LBM simulation functions to the list
void add_function(comp_func f, const char* name, int flops) {
    lbmFuncs.push_back(f);
    funcNames.push_back(name);
    funcFlops.push_back(flops);
    numFuncs++;
}

//Forward declaration of functions defined in lbm2d_baseline.c
extern "C" void register_functions();
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

    register_functions();
    if (numFuncs == 0) {
        cout << "No functions registered - exiting." << endl;
        return 0;
    }else{
        cout << endl << "[1/] Starting program. " << numFuncs << " function(s) registered." << endl;
    }

    /*-------------------------------------Benchmarking the Functions-------------------------------------------*/
    
    //iterate over various problem sizes, initialize F and cylinder, and run all functions. 
    //Nx = Ny
    cout << "[2/2] Testing correctness:" << endl;
    int sizes[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    int numSizes = 10;
    //array of pointers to the various F and cylinder arrays
    double* F[numSizes];
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
