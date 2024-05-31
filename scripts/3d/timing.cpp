#include <iostream>
#include <utility>
#include <vector>
#include <functional>
#include <fstream>
#include "tsc_x86.h"
#include "utils.h"

using namespace std;

#define CYCLES_REQUIRED 1e8
#define C_S 0.6
#define TAU 0.75
#define GAMMA_DOT 0.01

#define REP 10
#define TIME_STEPS 10

#define N_X 10
#define N_Y 15
#define N_Z 20
#define BOUNDARY_CONDITION_TYPE_DEFAULT 3
#define DIRECTIONS_TYPE_DEFAULT 27

//initialisation of BC is irrelevant for the testing and timing program

// due to current circumstances local variables are cumbersome to use ...
int DIRECTIONS_TYPE_g;
int NX_g;
int NY_g;
int NZ_g;
int BOUNDARY_CONDITION_g;


template <typename T>
class FuncEntry {
public:
    T func;
    std::function<void(T, LBMarrays*)> run_func;
    std::function<myInt64(T, int)> time_func;
    const char* funcName;
    std::function<struct ops(LBMarrays*)> calc_ops;
    FuncEntry(T p1, std::function<void(T, LBMarrays*)> p2, std::function<myInt64(T, int)> p3, const char* p4, std::function<struct ops(LBMarrays*)> p5) {
        func = p1;
        run_func = p2;
        time_func = p3;
        funcName = p4;
        calc_ops = std::move(p5);
    }
};

template <typename T>
double time_function(FuncEntry<T> f) {
    long num_runs = 1;
    double multiplier = 1;
    myInt64 time;
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do{
        num_runs = (long) ((double) num_runs * multiplier);
        auto cycles = (double) f.time_func(f.func, num_runs);
        multiplier = (CYCLES_REQUIRED) / cycles;
    } while (multiplier > 2);
    // Perform REP runs and store the results in a list
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        time = f.time_func(f.func, num_runs);
        double cycles = ((double)time) / (double) num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return total_cycles; // Returns the average cycles over the repetitions
}

int check_equality(LBMarrays* solver1, LBMarrays* solver2) {
    int nX = solver1->nX;
    int nY = solver1->nY;
    int nZ = solver1->nZ;
    int direction_size = solver1->direction_size;

    if(nX != solver2->nX || nY != solver2->nY || nZ != solver2->nZ || direction_size != solver2->direction_size) {
        cout << " dimensions of solvers do not match:" << endl;
        cout << "           expected x:" << nX << ", y: " << nY << ", z: " << nZ << ", directions " << direction_size << endl;
        cout << "           got      x:" << solver2->nX << ", y: " << solver2->nY << ", z: " << solver2->nZ << ", directions " << solver2->direction_size;
        return 0;
    }

    int eqPD = equal(solver1->particle_distributions, solver2->particle_distributions, nX * nY * nZ * direction_size);
    int eqPPD = equal(solver1->previous_particle_distributions, solver2->previous_particle_distributions, nX * nY * nZ * direction_size);
    int eqDF = equal(solver1->density_field, solver2->density_field, nX * nY * nZ);
    int eqVF = equalVelocity(solver1, solver2);

    if(!eqPD || !eqPPD || !eqDF || !eqVF) {
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
    }
    return eqPD + eqPPD + eqDF + eqVF;
}

// currently Testing Behaviour is not influenced by initialized boundary condition
struct LBMarrays* init_struct(int direction_size, int boundary_condition) {
    auto solver = (LBMarrays*) malloc(sizeof(LBMarrays));
    solver->nX = NX_g;
    solver->nY = NY_g;
    solver->nZ = NZ_g;

    solver->nXY = NX_g * NY_g;
    solver->nXYZ = NX_g * NY_g * NZ_g;

    solver->boundary_condition = boundary_condition; // 1=periodic, 2=couette, 3=lees_edwards
    solver->direction_size = direction_size; // One of 9, 15, 27

    solver->c_s = C_S;
    solver->tau = TAU;
    solver->gamma_dot = GAMMA_DOT;

    int box_flatten_length = solver->nXYZ;
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
    solver->velocity_fieldX = (double*) malloc(box_flatten_length * sizeof(double));
    solver->velocity_fieldY = (double*) malloc(box_flatten_length * sizeof(double));
    solver->velocity_fieldZ = (double*) malloc(box_flatten_length * sizeof(double));
    solver->previous_particle_distributions = (double*) malloc(distributions_flatten_length * sizeof(double));
    solver->particle_distributions = (double*) malloc(distributions_flatten_length * sizeof(double));
    solver->reverse_indexes = (int*) malloc(solver->direction_size * sizeof(int));

    initialize(solver->nX, solver->nY, solver->nZ, solver->direction_size, solver->density_field, solver->velocity_field, solver->previous_particle_distributions, solver->particle_distributions, solver->directions, solver->weights, solver->reverse_indexes);

    for(int i = 0; i < box_flatten_length; i++) {
        solver->velocity_fieldX[i] = solver->velocity_field[3 * i];
        solver->velocity_fieldY[i] = solver->velocity_field[3 * i + 1];
        solver->velocity_fieldZ[i] = solver->velocity_field[3 * i + 2];
    }

    return solver;
}


void free_struct(LBMarrays* solver) {
    free(solver->density_field);
    free(solver->velocity_field);
    free(solver->velocity_fieldX);
    free(solver->velocity_fieldY);
    free(solver->velocity_fieldZ);
    free(solver->previous_particle_distributions);
    free(solver->particle_distributions);
    free(solver->reverse_indexes);
    free(solver);
}

void run_func_struct(comp_func_struct f, LBMarrays* solver) {
    for(int j = 0; j < TIME_STEPS; j++) {
        f(solver);
    }
}

myInt64 time_func_struct(comp_func_struct f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays* solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_func_struct(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_func_struct_time(comp_func_struct_time f, LBMarrays* solver) {
    for(int j = 0; j < TIME_STEPS; j++) {
        f(solver, j);
    }
}

myInt64 time_func_struct_time(comp_func_struct_time f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays* solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_func_struct_time(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_collision_func_array(comp_collision_arrays f, LBMarrays* solver) {
    for (int j = 0; j < TIME_STEPS; j++) {
        f(
                solver->nX, solver->nY, solver->nZ, solver->direction_size, solver->tau, solver->c_s,
                solver->density_field,
                solver->velocity_field,
                solver->previous_particle_distributions,
                solver->particle_distributions,
                solver->directions,
                solver->weights
        );
    }
}

myInt64 time_collision_func_array(comp_collision_arrays f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_collision_func_array(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_momentum_func_array(comp_momentum_arrays f, LBMarrays* solver) {
    for (int j = 0; j < TIME_STEPS; j++) {
        f(
                solver->nX, solver->nY, solver->nZ, solver->direction_size,
                solver->density_field,
                solver->velocity_field,
                solver->particle_distributions,
                solver->directions
        );
    }
}

myInt64 time_momentum_func_array(comp_momentum_arrays f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_momentum_func_array(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_stream_periodic_func_array(comp_stream_periodic_arrays f, LBMarrays* solver) {
    for (int j = 0; j < TIME_STEPS; j++) {
        f(
                solver->nX, solver->nY, solver->nZ, solver->direction_size,
                solver->previous_particle_distributions,
                solver->particle_distributions,
                solver->directions
        );
    }
}

myInt64 time_stream_periodic_func_array(comp_stream_periodic_arrays f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_stream_periodic_func_array(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_stream_couette_func_array(comp_stream_couette_arrays f, LBMarrays* solver) {
    for (int j = 0; j < TIME_STEPS; j++) {
        f(
                solver->nX, solver->nY, solver->nZ, solver->direction_size,
                solver->c_s,
                solver->previous_particle_distributions,
                solver->particle_distributions,
                solver->directions,
                solver->weights,
                solver->reverse_indexes
        );
    }
}

myInt64 time_stream_couette_func_array(comp_stream_couette_arrays f, long num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_stream_couette_func_array(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_stream_lees_edwards_func_array(comp_stream_lees_edwards_arrays f, LBMarrays* solver) {
    for (int j = 0; j < TIME_STEPS; j++) {
        f(
                solver->nX, solver->nY, solver->nZ, solver->direction_size, j, solver->gamma_dot,
                solver->c_s,
                solver->density_field,
                solver->velocity_field,
                solver->previous_particle_distributions,
                solver->particle_distributions,
                solver->directions,
                solver->weights
        );
    }
}

myInt64 time_stream_lees_edwards_func_array(comp_stream_lees_edwards_arrays f, int num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_stream_lees_edwards_func_array(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

void run_func_array(comp_func_arrays f, LBMarrays* solver) {
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
}

myInt64 time_func_array(comp_func_arrays f, int num_runs) {
    myInt64 start;
    myInt64 total = 0;
    for(int i = 0; i < num_runs; i++) {
        LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        start = start_tsc();
        run_func_array(f, solver);
        total += stop_tsc(start);
        free_struct(solver);
    }
    return total;
}

template <typename T, typename U>
void step(int num, int max, std::ofstream& fos, const char* name, FuncEntry<T> baseline, vector<FuncEntry<T>> structFuncs, vector<FuncEntry<U>> arrayFuncs) {
    struct LBMarrays *example = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);

    size_t numFuncs = structFuncs.size() + arrayFuncs.size();
    if (numFuncs == 0) {
        cout << endl << "[" << num << "/" << max << "] " << name << ": No functions registered, skipping..." << endl;
    } else {
        cout << endl << "[" << num << "/" << max << "] " << name << ": " << numFuncs << " function(s) registered." << endl;
        cout << "    [1/2] Testing correctness:" << endl;
        struct LBMarrays *baseline_solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
        baseline.run_func(baseline.func, baseline_solver);
        int eq = 0;
        for (int i = 0; i < structFuncs.size(); i++) {
            // Run the function to be tested
            cout << "       Testing [" << i + 1 << "/" << numFuncs << "]: Function \"" << structFuncs[i].funcName
                 << "\"" << endl;
            struct LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
            structFuncs[i].run_func(structFuncs[i].func, solver);
            eq += check_equality(baseline_solver, solver);
            free_struct(solver);
        }
        for (int i = 0; i < arrayFuncs.size(); i++) {
            // Run the function to be tested
            cout << "       Testing [" << structFuncs.size() + i + 1 << "/" << numFuncs << "]: Function \""
                 << arrayFuncs[i].funcName << "\"" << endl;
            struct LBMarrays *solver = init_struct(DIRECTIONS_TYPE_g, BOUNDARY_CONDITION_g);
            arrayFuncs[i].run_func(arrayFuncs[i].func, solver);
            eq += check_equality(baseline_solver, solver);
            free_struct(solver);
        }
        free_struct(baseline_solver);
        if (eq < numFuncs * 4) {
            cout << eq;
            cout << "       Testing FAILED: Some simulations compute incorrect results." << endl;
            return;
        } else {
            cout << "    [2/2] Timing:" << endl;
            double cycles = time_function(baseline);
            struct ops baselineOps = baseline.calc_ops(example);
            fos << baselineOps.iops   << "," <<
                   baselineOps.flops  << "," <<
                   baselineOps.bytes_read  << "," <<
                   baselineOps.bytes_write  << "," <<
                   name << "," <<
                   baseline.funcName               << "," <<
                   cycles                          << "," <<
                   example->direction_size         << "," <<
                   NX_g                             << "," <<
                   NY_g                             << "," <<
                   NZ_g                             << "," <<
                   TIME_STEPS                      << "," <<
                   std::endl;
            cout << "       Timing [1/" << numFuncs + 1 << "]: Baseline on average takes " << cycles << " cycles for "
                 << TIME_STEPS << " time step." << endl;
            for (int i = 0; i < structFuncs.size(); i++) {
                cycles = time_function(structFuncs[i]);
                cout << "       Timing [" << i + 2 << "/" << numFuncs + 1 << "]: Function \"" << structFuncs[i].funcName
                     << "\" on average takes " << cycles << " cycles for " << TIME_STEPS << " time step." << endl;
                struct ops ops = structFuncs[i].calc_ops(example);
                fos << ops.iops   << "," <<
                       ops.flops  << "," <<
                       ops.bytes_read  << "," <<
                       ops.bytes_write  << "," <<
                       name << "," <<
                       structFuncs[i].funcName << "," <<
                       cycles                  << "," <<
                       example->direction_size << "," <<
                       NX_g                     << "," <<
                       NY_g                     << "," <<
                       NZ_g                     << "," <<
                       TIME_STEPS              << "," <<
                       std::endl;
            }
            for (int i = 0; i < arrayFuncs.size(); i++) {
                cycles = time_function(arrayFuncs[i]);
                cout << "       Timing [" << structFuncs.size() + i + 2 << "/" << numFuncs + 1 << "]: Function \""
                     << arrayFuncs[i].funcName << "\" on average takes " << cycles << " cycles for " << TIME_STEPS
                     << " time step." << endl;

                struct ops ops = arrayFuncs[i].calc_ops(example);
                fos << ops.iops   << "," <<
                       ops.flops  << "," <<
                       ops.bytes_read  << "," <<
                       ops.bytes_write  << "," <<
                       name << "," <<
                       arrayFuncs[i].funcName  << "," <<
                       cycles                  << "," <<
                       example->direction_size << "," <<
                       NX_g                    << "," <<
                       NY_g                    << "," <<
                       NZ_g                    << "," <<
                       TIME_STEPS              << "," <<
                       std::endl;
            }
        }
        cout << "       Timing DONE." << endl;
        free_struct(example);
    }
}

// Global vars, used to keep track of LBM simulation functions
vector<FuncEntry<comp_func_struct>> momentumFuncsStruct;
vector<FuncEntry<comp_func_struct>> collisionFuncsStruct;
vector<FuncEntry<comp_func_struct_time>> streamPeriodicFuncsStruct;
vector<FuncEntry<comp_func_struct>> streamCouetteFuncsStruct;
vector<FuncEntry<comp_func_struct_time>> streamLeesEdwardsFuncsStruct;
vector<FuncEntry<comp_func_struct_time>> lbmFuncsStruct;

vector<FuncEntry<comp_momentum_arrays>> momentumFuncsArrays;
vector<FuncEntry<comp_collision_arrays>> collisionFuncsArrays;
vector<FuncEntry<comp_stream_periodic_arrays>> streamPeriodicFuncsArrays;
vector<FuncEntry<comp_stream_couette_arrays>> streamCouetteFuncsArrays;
vector<FuncEntry<comp_stream_lees_edwards_arrays>> streamLeesEdwardsFuncsArrays;
vector<FuncEntry<comp_func_arrays>> lbmFuncsArrays;

// Function to add LBM simulation functions to the lists
void add_momentum_struct_func(comp_func_struct f, calc_flops calc_ops, const char* name) {
    momentumFuncsStruct.emplace_back(f, &run_func_struct, &time_func_struct, name, calc_ops);
}
void add_collision_struct_func(comp_func_struct f, calc_flops calc_ops, const char* name) {
    collisionFuncsStruct.emplace_back(f, &run_func_struct, &time_func_struct, name, calc_ops);
}
void add_stream_periodic_struct_func(comp_func_struct_time f, calc_flops calc_ops, const char* name) {
    streamPeriodicFuncsStruct.emplace_back(f, &run_func_struct_time, &time_func_struct_time, name, calc_ops);
}
void add_stream_couette_struct_func(comp_func_struct f, calc_flops calc_ops, const char* name) {
    streamCouetteFuncsStruct.emplace_back(f, &run_func_struct, &time_func_struct, name, calc_ops);
}
void add_stream_lees_edwards_struct_func(comp_func_struct_time f, calc_flops calc_ops, const char* name) {
    streamLeesEdwardsFuncsStruct.emplace_back(f, &run_func_struct_time, &time_func_struct_time, name, calc_ops);
}
void add_lbm_struct_func(comp_func_struct_time f, calc_flops calc_ops, const char* name) {
    lbmFuncsStruct.emplace_back(f, &run_func_struct_time, &time_func_struct_time, name, calc_ops);
}

void add_momentum_array_func(comp_momentum_arrays f, calc_flops calc_ops, const char* name) {
    momentumFuncsArrays.emplace_back(f, &run_momentum_func_array, &time_momentum_func_array, name, calc_ops);
}
void add_collision_array_func(comp_collision_arrays f, calc_flops calc_ops, const char* name) {
    collisionFuncsArrays.emplace_back(f, &run_collision_func_array, &time_collision_func_array, name, calc_ops);
}
void add_stream_periodic_array_func(comp_stream_periodic_arrays f, calc_flops calc_ops, const char* name) {
    streamPeriodicFuncsArrays.emplace_back(f, &run_stream_periodic_func_array, &time_stream_periodic_func_array, name, calc_ops);
}
void add_stream_couette_array_func(comp_stream_couette_arrays f, calc_flops calc_ops, const char* name) {
    streamCouetteFuncsArrays.emplace_back(f, &run_stream_couette_func_array, &time_stream_couette_func_array, name, calc_ops);
}
void add_stream_lees_edwards_array_func(comp_stream_lees_edwards_arrays f, calc_flops calc_ops, const char* name) {
    streamLeesEdwardsFuncsArrays.emplace_back(f, &run_stream_lees_edwards_func_array, &time_stream_lees_edwards_func_array, name, calc_ops);
}
void add_lbm_array_func(comp_func_arrays f, calc_flops calc_ops, const char* name) {
    lbmFuncsArrays.emplace_back(f, &run_func_array, &time_func_array, name, calc_ops);
}

int main(int argc, char* argv[]) {
    bool test_collision(true);
    bool test_momentum(true);
    int test_stream(4);
    bool test_LBM(false);
    bool reset_datafile(false);
    NX_g = N_X;
    NY_g = N_Y;
    NZ_g = N_Z;
    DIRECTIONS_TYPE_g = DIRECTIONS_TYPE_DEFAULT;
    switch(argc-1){
        // case 5:
        //     BOUNDARY_CONDITION_g = std::stoi(argv[3]);
        case 9:
            reset_datafile = (bool)(std::stoi(argv[9]));
        case 8:
            test_LBM = (bool)(std::stoi(argv[8]));
        case 7:
            test_collision = (bool)(std::stoi(argv[7]));
        case 6:
            test_momentum = (bool)(std::stoi(argv[6]));
        case 5:
            test_stream = std::stoi(argv[5]); //0 none; 1 periodic 2 couette 3 lees edwards 4 all
        case 4:
            NZ_g = std::stoi(argv[4]);
        case 3:
            NY_g = std::stoi(argv[3]);
        case 2:
            NX_g = std::stoi(argv[2]);
        case 1:
            DIRECTIONS_TYPE_g = std::stoi(argv[1]);
        case 0: break;
        default: std::cout << "too many arguments, default initialisation will be taken" << std::endl;
    }

    // register_functions();
    register_momentum_functions();
    register_collision_functions();
    register_stream_periodic_functions();
    register_stream_couette_functions();
    register_stream_lees_edwards_functions();
    register_lbm_function();

    std::string filename = "TimingData.csv";
    std::ofstream fos;

    if(reset_datafile) {
        fos.open(filename, std::ofstream::out);
        fos << "iops,flops,bytes_read,bytes_write,step,function,cycles,DIRECTION_SIZE,NX,NY,NZ,TIMESTEPS,bytes" << std::endl;
    } else {
        fos.open(filename, std::ofstream::out | std::ofstream::app);
    }

    if(test_momentum) {
        auto momentumBaseline = FuncEntry<comp_func_struct>(&momentum_baseline, &run_func_struct, &time_func_struct, "Momentum Baseline", &momentum_baseline_flops);
        step(1, 8, fos, "Momentum", momentumBaseline, momentumFuncsStruct, momentumFuncsArrays);
    }
    if(test_collision) {
        auto collisionBaseline = FuncEntry<comp_func_struct>(&collision_baseline, &run_func_struct, &time_func_struct, "Collision Baseline", &collision_baseline_flops);
        step(2, 8, fos, "Collision", collisionBaseline, collisionFuncsStruct, collisionFuncsArrays);
    }
    if(test_stream == 1 || test_stream == 4) {
        auto streamPeriodicBaseline = FuncEntry<comp_func_struct_time>(&stream_periodic_baseline, &run_func_struct_time, &time_func_struct_time, "Stream Periodic Baseline", &stream_periodic_baseline_flops);
        step(3, 8, fos, "Stream Periodic", streamPeriodicBaseline, streamPeriodicFuncsStruct, streamPeriodicFuncsArrays);
    }
    if(test_stream == 2 || test_stream == 4) {
        auto streamCouetteBaseline = FuncEntry<comp_func_struct>(&stream_couette_baseline, &run_func_struct, &time_func_struct, "Stream Couette Baseline", &stream_couette_baseline_flops);
        step(4, 8, fos, "Stream Couette", streamCouetteBaseline, streamCouetteFuncsStruct, streamCouetteFuncsArrays);
    }
    if(test_stream == 3 || test_stream == 4) {
        auto streamLeesEdwardsBaseline = FuncEntry<comp_func_struct_time>(&stream_lees_edwards_baseline, &run_func_struct_time, &time_func_struct_time, "Stream Lees Edwards Baseline", &stream_lees_edwards_baseline_flops);
        step(5, 8, fos, "Stream Less Edwards", streamLeesEdwardsBaseline, streamLeesEdwardsFuncsStruct, streamLeesEdwardsFuncsArrays);
    }
    if(test_LBM){
        auto lbmBaseline = FuncEntry<comp_func_struct_time>(&perform_timestep_baseline, &run_func_struct_time, &time_func_struct_time, "LBM Baseline", &perform_timestep_baseline_flops);
        BOUNDARY_CONDITION_g = 1;
        step(6, 8, fos, "LBM - Periodic Boundary Condition", lbmBaseline, lbmFuncsStruct, lbmFuncsArrays);
        BOUNDARY_CONDITION_g = 2;
        step(7, 8, fos, "LBM - Couette Boundary Condition", lbmBaseline, lbmFuncsStruct, lbmFuncsArrays);
        BOUNDARY_CONDITION_g = 3;
        step(8, 8, fos, "LBM - Lees Edwards Boundary Condition", lbmBaseline, lbmFuncsStruct, lbmFuncsArrays);
    }

    cout << "Program finished." << endl << endl;

    fos.close();
    return 0;
}
