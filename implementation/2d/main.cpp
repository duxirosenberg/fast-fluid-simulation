#include <iostream>
#include <vector>
#include <string>
#include <list>
#include "tsc_x86.h"
#include "utils.h"

using namespace std;

#define CYCLES_REQUIRED 1e8
#define REP 30
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

//Forward declaration of functions defined in lbm_2d.c
extern "C" void register_functions();
extern "C" void lbm_baseline(double *F, int *cylinder, int Nx, int Ny, int n_steps);


double time_function(comp_func f, int flops) {
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;

    // Example problem size and data
    //ToDo Douglas: add a loop to test different problem sizes
    const int Nx = 10, Ny = 10;
    double* F = (double*)malloc(Nx*Ny*NL * sizeof(double));
    int* cylinder = (int*)malloc(Nx*Ny * sizeof(int));

    // Build initial conditions for F and cylinder
    initialize_F(F, Nx, Ny);
    initialize_cylinder(cylinder, Nx, Ny);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.  
    do{
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            //perform Nt steps of Simulation on initial conditions F_init, store the result in F
            f(F, cylinder, Nx, Ny, Nt );
        }
        end = stop_tsc(start);

        cycles = (double) end;
        multiplier = (CYCLES_REQUIRED) / cycles;
    } while (multiplier > 2);

    list<double> cyclesList;

    //reinitialize F
    initialize_F(F, Nx, Ny);
    // Perform REP runs and store the results in a list
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            //perform Nt steps of Simulation on initial conditions F_init, store the result in F
            f(F, cylinder, Nx, Ny, Nt);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;
        cyclesList.push_back(cycles);
    }
    total_cycles /= REP;

    cyclesList.sort();

    free(F);
    free(cylinder);

    return total_cycles; // Returns the average cycles over the repetitions
}

int main(int argc, char **argv) {

    register_functions();
    if (numFuncs == 0) {
        cout << "No functions registered - exiting." << endl;
        return 0;
    }else{
        cout << endl << "[1/3] Starting program. " << numFuncs << " function(s) registered." << endl;
    }


    /*-------------------------------------Validating the Functions-------------------------------------------*/
    cout << "[2/3] Testing correctness:" << endl;
    int Nx = 400, Ny = 100, n_steps = 4000;
    double* F = (double*)calloc(Nx*Ny*NL,sizeof(double));
    double* F_baseline = (double*)calloc(Nx*Ny*NL , sizeof(double));
    int* cylinder = (int*)calloc(Nx*Ny , sizeof(int));

    // initialize the cylinder
    initialize_cylinder(cylinder, Nx, Ny);
   
    // Run the baseline to compare the other functions to
    cout << "       Testing [1/" << numFuncs+1 << "]: Currently Calculating Baseline " << endl;
    initialize_F(F_baseline, Nx, Ny);
    lbm_baseline(F_baseline, cylinder, Nx, Ny, n_steps);

    for (int i = 0; i < numFuncs; i++) {
        //reset F
        F = (double*)calloc(Nx*Ny*NL,sizeof(double));

        // Run the function to be tested
        cout << "       Testing [" << i+2 << "/" << numFuncs+1 << "]: Function \"" << funcNames[i] << "\"";
        initialize_F(F, Nx, Ny);
        lbmFuncs[i](F, cylinder, Nx, Ny, n_steps);
        
        // Check the result
        int eq = equal(F_baseline, F, Nx*Ny*NL);
        if (eq) {
            cout << " computes CORRECT result. "<< endl;

        } else {
            cout << " computes INCORRECT result. "<< endl;
            exit(1);
        }
    }
    cout << "       Testing DONE: All simulations compute correct results." << endl;


    /*-------------------------------------Timing the Functions-------------------------------------------*/
    cout << "[3/3] Timing the functions." << endl;
    double cycles = time_function(lbm_baseline, 0);
    cout << "       Timing [1/" << numFuncs+1 << "]: Baseline on average takes " << cycles << " cycles per timestep." << endl;
    for (int i = 0; i < numFuncs; i++) {
        cycles = time_function(lbmFuncs[i], funcFlops[i]);
        cout << "       Timing [" << i+2 << "/" << numFuncs+1 << "]: Function \"" << funcNames[i] << "\" on average takes " << cycles << " cycles per timestep." <<endl;
    }
    cout << "       Timing DONE." << endl;
    cout << "Program finished." << endl<<endl;

    return 0;
}
