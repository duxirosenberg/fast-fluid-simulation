#include <iostream>
#include <vector>
#include <string>
#include <list>
#include "tsc_x86.h"
#include "utils2d.h"

using namespace std;

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


int main(int argc, char **argv) {

    register_functions();//initializes global vars lbmFuncs, funcNames, funcFlops, numFuncs
    if (numFuncs == 0) {
        cout << "No functions registered - exiting." << endl;
        return 0;
    }else{
        cout << endl << "[1/2] Starting program. " << numFuncs << " function(s) registered." << endl;
    }


    /*-------------------------------------Validating the Functions-------------------------------------------*/
    cout << "[2/2] Testing correctness:" << endl;
    int Nx = 400, Ny = 100, n_steps = 10;
    double* F = (double*)calloc(Nx*Ny*NL,sizeof(double));
    double* F_baseline = (double*)calloc(Nx*Ny*NL , sizeof(double));
    int* cylinder = (int*)calloc(Nx*Ny , sizeof(int));

    // initialize the cylinder
    initialize_cylinder(cylinder, Nx, Ny);
   
    // Run the baseline to compare the other functions to
    cout << "       Testing [0/" << numFuncs+1 << "]: Initializing Conditions" << endl;
    initialize_F(F_baseline, Nx, Ny);
    cout << "       Testing [1/" << numFuncs+1 << "]: Running Baseline Simulation" << endl;
    lbm_baseline(F_baseline, cylinder, Nx, Ny, n_steps);
    free(F);

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
        //free
        free(F);
    }
    //free mem
    free(F_baseline);
    free(cylinder);
    cout << "       Testing DONE: All simulations compute correct results." << endl;

    return 0;
}
