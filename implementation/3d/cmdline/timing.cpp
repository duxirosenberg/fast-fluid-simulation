#include <iostream>
#include <vector>
#include <string>
#include "tsc_x86.h"
#include "utils.h"

using namespace std;

#define CYCLES_REQUIRED 1e9
#define REP 10
#define TIME_STEPS 20
#define N_X 10
#define N_Y 15
#define N_Z 20
#define C_S 0.6
#define TAU 0.75
#define GAMMA_DOT 0.01
#define TIME_STEPS 20
#define EPS 1e-3

// Global vars, used to keep track of LBM simulation functions
vector<comp_func_struct> lbmFuncsStruct;
vector<comp_func_arrays> lbmFuncsArrays;
vector<string> funcNamesArrays;
vector<string> funcNamesStruct;
vector<int> funcFlopsStruct;
vector<int> funcFlopsArrays;

// Function to add LBM simulation functions to the list
void add_array_func(comp_func_arrays f, const char* name, int flops) {
    lbmFuncsArrays.push_back(f);
    funcNamesArrays.emplace_back(name);
    funcFlopsArrays.push_back(flops);
}

void add_struct_func(comp_func_struct f, const char* name, int flops) {
    lbmFuncsStruct.push_back(f);
    funcNamesStruct.emplace_back(name);
    funcFlopsStruct.push_back(flops);
}

struct LBMarrays* init_struct(int direction_size, int boundary_condition) {
    auto solver = (LBMarrays*) malloc(sizeof(LBMarrays));
    solver->nX = N_X;
    solver->nY = N_Y;
    solver->nZ = N_Z;
    solver->direction_size = direction_size; // One of 9, 15, 27
    solver->c_s = C_S;
    solver->tau = TAU;
    solver->gamma_dot = GAMMA_DOT;
    solver->boundary_condition = boundary_condition; // 1=periodic, 2=couette, 3=lees_edwards

    int box_flatten_length = solver->nX * solver->nY * solver->nZ;
    int distributions_flatten_length = box_flatten_length * solver->direction_size;
    if(solver->direction_size == 15) {
        solver->directions = &D3Q15_DIRECTIONS[0];
        solver->weights = &D3Q15_WEIGHTS[0];
    } else if(solver->direction_size == 27) {
        solver->directions = &D3Q27_DIRECTIONS[0];
        solver->weights = &D3Q27_WEIGHTS[0];
    } else if(solver->direction_size == 9) {
        // 2D case
        solver->directions = &D2Q9_DIRECTIONS[0];
        solver->weights = &D2Q9_WEIGHTS[0];
    }
    solver->density_field = (double*) malloc(box_flatten_length * sizeof(double));
    solver->velocity_field = (double*) malloc(3 * box_flatten_length * sizeof(double));
    solver->previous_particle_distributions = (double*) malloc(distributions_flatten_length * sizeof(double));
    solver->particle_distributions = (double*) malloc(distributions_flatten_length * sizeof(double));
    solver->reverse_indexes = (int*) malloc(solver->direction_size * sizeof(int));

    initialize(solver->nX, solver->nY, solver->nZ, solver->direction_size, solver->density_field, solver->velocity_field, solver->previous_particle_distributions, solver->particle_distributions, solver->directions, solver->weights, solver->reverse_indexes);
    return solver;
}


void free_struct(LBMarrays* solver) {
    free(solver->density_field);
    free(solver->velocity_field);
    free(solver->previous_particle_distributions);
    free(solver->particle_distributions);
    free(solver->reverse_indexes);
    free(solver);
}

myInt64 time_func_struct(comp_func_struct f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays* solver = init_struct(27, 3);
        start = start_tsc();
        for(int j = 0; j < TIME_STEPS; j++) {
            f(solver, j);
        }
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}


myInt64 time_func_array(comp_func_arrays f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(27, 3);
        start = start_tsc();
        for (int j = 0; j < TIME_STEPS; j++) {
            f(
                    solver->nX, solver->nY, solver->nZ, solver->direction_size, j, solver->tau, solver->gamma_dot,
                    solver->c_s, solver->boundary_condition,
                    solver->density_field,
                    solver->velocity_field,
                    solver->previous_particle_distributions,
                    solver->particle_distributions,
                    solver->directions,
                    solver->weights,
                    solver->reverse_indexes
            );
        }
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

double time_function_struct(comp_func_struct f, int flops) {
    long num_runs = 20;
    double multiplier = 1;
    myInt64 time;
    // Build initial conditions for F and cylinder
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do{
        num_runs = (long) ((double) num_runs * multiplier);
        // perform Nt steps of Simulation on initial conditions F_init, store the result in F
        auto cycles = (double) time_func_struct(f, num_runs);
        multiplier = (CYCLES_REQUIRED) / cycles;
    } while (multiplier > 2);
    // Perform REP runs and store the results in a list
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        time = time_func_struct(f, num_runs);
        double cycles = ((double)time) / (double) num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return total_cycles; // Returns the average cycles over the repetitions
}

double time_function_array(comp_func_arrays f, int flops) {
    long num_runs = 20;
    double multiplier = 1;
    myInt64 time;
    // Build initial conditions for F and cylinder
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do{
        num_runs = (long) ((double) num_runs * multiplier);
        // perform Nt steps of Simulation on initial conditions F_init, store the result in F
        auto cycles = (double) time_func_array(f, num_runs);
        multiplier = (CYCLES_REQUIRED) / cycles;
    } while (multiplier > 2);
    // Perform REP runs and store the results in a list
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        time = time_func_array(f, num_runs);
        double cycles = ((double)time) / (double) num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return total_cycles; // Returns the average cycles over the repetitions
}


int check_equality(LBMarrays* solver1, LBMarrays* solver2) {
    int nX = solver1->nX;
    int nY = solver1->nX;
    int nZ = solver1->nZ;
    int direction_size = solver1->direction_size;

    if(nX != solver2->nX || nY != solver2->nX || nZ != solver2->nZ || direction_size != solver2->direction_size) {
        cout << " dimensions of solvers do not match:" << endl;
        cout << "           expected x:" << nX << ", y: " << nY << ", z: " << nZ << ", directions " << direction_size << endl;
        cout << "           got      x:" << solver2->nX << ", y: " << solver2->nY << ", z: " << solver2->nZ << ", directions " << solver2->direction_size;
        return 0;
    }

    int eqPD = equal(solver1->particle_distributions, solver2->particle_distributions, nX * nY * nZ * direction_size);
    int eqPPD = equal(solver1->previous_particle_distributions, solver2->previous_particle_distributions, nX * nY * nZ * direction_size);
    int eqDF = equal(solver1->density_field, solver2->density_field, nX * nY * nZ);
    int eqVF = equal(solver1->velocity_field, solver2->velocity_field, nX * nY * nZ);

    if (eqPD) {
        cout << "           particle_distributions: CORRECT" << endl;
    } else {
        cout << "           particle_distributions: INCORRECT" << endl;
    }
    if (eqPPD) {
        cout << "           previous_particle_distributions: CORRECT" << endl;
    } else {
        cout << "           previous_particle_distributions: INCORRECT" << endl;
    }
    if (eqDF) {
        cout << "           density_field: CORRECT" << endl;
    } else {
        cout << "           density_field: INCORRECT" << endl;
    }
    if (eqVF) {
        cout << "           velocity_field: CORRECT" << endl;
    } else {
        cout << "           velocity_field: INCORRECT" << endl;
    }
    return eqPD + eqPPD + eqDF + eqVF;
}

int main(int argc, char **argv) {
    register_functions();

    size_t numFuncs = lbmFuncsStruct.size() + lbmFuncsArrays.size() + 1;
    if (numFuncs == 0) {
        cout << "No functions registered - exiting." << endl;
        return 0;
    }else{
        cout << endl << "[1/3] Starting program. " << numFuncs << " function(s) registered." << endl;
    }


    /*-------------------------------------Validating the Functions-------------------------------------------*/
    cout << "[2/3] Testing correctness:" << endl;

    // initialize
    struct LBMarrays* baseline_solver = init_struct(27, 3);
    // Run the baseline to compare the other functions to
    cout << "       Testing [1/" << numFuncs << "]: Currently Calculating Baseline " << endl;
    for(int i = 0; i < TIME_STEPS; i++) {
        perform_timestep_baseline(baseline_solver, i);
    }

    int eq = 0;
    for (int i = 0; i < lbmFuncsStruct.size(); i++) {
        // Run the function to be tested
        cout << "       Testing [" << i + 2 << "/" << numFuncs << "]: Function \"" << funcNamesStruct[i] << "\"" << endl;
        struct LBMarrays* solver = init_struct(27, 3);
        for(int j = 0; j < TIME_STEPS; j++) {
            lbmFuncsStruct[i](solver, j);
        }
        eq += check_equality(baseline_solver, solver);
        free_struct(solver);
    }

    for (int i = 0; i < lbmFuncsArrays.size(); i++) {
        // Run the function to be tested
        cout << "       Testing [" << lbmFuncsArrays.size() + i + 2 << "/" << numFuncs << "]: Function \"" << funcNamesArrays[i] << "\"" << endl;
        struct LBMarrays* solver = init_struct(27, 3);
        for(int j = 0; j < TIME_STEPS; j++) {
            lbmFuncsArrays[i](
                    solver->nX, solver->nY, solver->nZ, solver->direction_size, j, solver->tau, solver->gamma_dot, solver->c_s, solver->boundary_condition,
                    solver->density_field,
                    solver->velocity_field,
                    solver->previous_particle_distributions,
                    solver->particle_distributions,
                    solver->directions,
                    solver->weights,
                    solver->reverse_indexes
            );
        }
        eq += check_equality(baseline_solver, solver);
        free_struct(solver);
    }
    free_struct(baseline_solver);
    if(eq < (numFuncs - 1) * 4) {
        cout << eq;
        cout << "       Testing FAILED: Some simulations compute incorrect results." << endl;
        exit(1);
    }
    cout << "       Testing DONE: All simulations compute correct results." << endl;


    /*-------------------------------------Timing the Functions-------------------------------------------*/
    cout << "[3/3] Timing the functions." << endl;
    double cycles = time_function_struct(&perform_timestep_baseline, 10);
    cout << "       Timing [1/" << numFuncs << "]: Baseline on average takes " << cycles << " cycles for " << TIME_STEPS << " timestep." << endl;
    for (int i = 0; i < lbmFuncsStruct.size(); i++) {
        cycles = time_function_struct(lbmFuncsStruct[i], funcFlopsStruct[i]);
        cout << "       Timing [" << i+2 << "/" << numFuncs << "]: Function \"" << funcNamesStruct[i] << "\" on average takes " << cycles << " cycles for " << TIME_STEPS << " timestep." << endl;
    }
    for (int i = 0; i < lbmFuncsArrays.size(); i++) {
        cycles = time_function_array(lbmFuncsArrays[i], funcFlopsArrays[i]);
        cout << "       Timing [" << lbmFuncsStruct.size() + i + 2 << "/" << numFuncs << "]: Function \"" << funcNamesArrays[i] << "\" on average takes " << cycles << " cycles for " << TIME_STEPS << " timestep." << endl;
    }
    cout << "       Timing DONE." << endl;
    cout << "Program finished." << endl << endl;

    return 0;
}
