#include <iostream>
#include <vector>
#include <string>
#include <list>

#include "tsc_x86.h"
#include "utils2d.h"
#include "lbm2d_struct.h"

using namespace std;

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
void add_drift_func(f_structs f, const char* name, int flops) {
    Funcs_drift.push_back(f);
    funcNames_drift.push_back(name);
    funcFlops_drift.push_back(flops);
    numFuncs_drift++;
}

// Function to add new reflective boundaries optimization
void add_reflective_boundaries_func(f_structs f, const char* name, int flops) {
    Funcs_set_reflective_boundaries.push_back(f);
    funcNames_set_reflective_boundaries.push_back(name);
    funcFlops_set_reflective_boundaries.push_back(flops);
    numFuncs_set_reflective_boundaries++;
}

// Function to add new collision optimization
void add_collision_func(f_structs f, const char* name, int flops) {
    Funcs_collision.push_back(f);
    funcNames_collision.push_back(name);
    funcFlops_collision.push_back(flops);
    numFuncs_collision++;
}

// Function to add new apply boundary optimization
void add_apply_boundary_func(f_structs f, const char* name, int flops) {
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
extern "C" struct LBM* lbm_initialize(int Nx, int Ny, const char* filename);
extern "C" void lbm_finalize(struct LBM* S);



int main(int argc, char **argv) {

    register_drift_functions();
    register_set_reflective_boundaries_functions();
    register_collision_functions();
    register_apply_boundary_functions();


    cout << endl << "[1/] Starting program. " << numFuncs_drift + numFuncs_set_reflective_boundaries + numFuncs_collision + numFuncs_apply_bndry << " function(s) registered." << endl;
    cout << "   Drift functions: " << numFuncs_drift << endl;
    cout << "   Reflective boundaries functions: " << numFuncs_set_reflective_boundaries << endl;
    cout << "   Collision functions: " << numFuncs_collision << endl;
    cout << "   Apply boundary functions: " << numFuncs_apply_bndry << endl;



    /*-------------------------------------Validating the Functions-------------------------------------------*/
    cout << "[2/] Testing correctness:" << endl;
    int Nx = 400, Ny = 100, n_steps = 10;

    //array of structs for each function
    struct LBM* test_drift_funcs[numFuncs_drift];
    struct LBM* test_set_reflective_boundaries_funcs[numFuncs_set_reflective_boundaries];
    struct LBM* test_collision_funcs[numFuncs_collision];
    struct LBM* test_apply_bndry_funcs[numFuncs_apply_bndry];

    //initialize all structs
    for(int i = 0; i<numFuncs_drift; i++){
        test_drift_funcs[i] = lbm_initialize(Nx, Ny, "");
    }
    for(int i = 0; i<numFuncs_set_reflective_boundaries; i++){
        test_set_reflective_boundaries_funcs[i] = lbm_initialize(Nx, Ny, "");
    }
    for(int i = 0; i<numFuncs_collision; i++){
        test_collision_funcs[i] = lbm_initialize(Nx, Ny, "");
    }
    for(int i = 0; i<numFuncs_apply_bndry; i++){
        test_apply_bndry_funcs[i] = lbm_initialize(Nx, Ny, "");
    }

    //apply functions to each struct, compare against baseline (func[0]) for each
    for(int i = 0; i<numFuncs_drift; i++){
        Funcs_drift[i](test_drift_funcs[i]);
        if(equal(test_drift_funcs[i]->F, test_drift_funcs[0]->F, Nx*Ny*NL)){
            cout << "   Function \"" << funcNames_drift[i] << "\" computes CORRECT result." << endl;
        }else{
            cout << "   Function \"" << funcNames_drift[i] << "\" computes INCORRECT result." << endl;
            exit(1);
        }
    }
    for(int i = 0; i<numFuncs_set_reflective_boundaries; i++){
        Funcs_set_reflective_boundaries[i](test_set_reflective_boundaries_funcs[i]);
        if(equal(test_set_reflective_boundaries_funcs[i]->F, test_set_reflective_boundaries_funcs[0]->F, Nx*Ny*NL)){
            cout << "   Function \"" << funcNames_set_reflective_boundaries[i] << "\" computes CORRECT result." << endl;
        }else{
            cout << "   Function \"" << funcNames_set_reflective_boundaries[i] << "\" computes INCORRECT result." << endl;
            exit(1);
        }
    }
    for(int i = 0; i<numFuncs_collision; i++){
        Funcs_collision[i](test_collision_funcs[i]);
        if(equal(test_collision_funcs[i]->F, test_collision_funcs[0]->F, Nx*Ny*NL)){
            cout << "   Function \"" << funcNames_collision[i] << "\" computes CORRECT result." << endl;
        }else{
            cout << "   Function \"" << funcNames_collision[i] << "\" computes INCORRECT result." << endl;
            exit(1);
        }
    }
    for(int i = 0; i<numFuncs_apply_bndry; i++){
        Funcs_apply_bndry[i](test_apply_bndry_funcs[i]);
        if(equal(test_apply_bndry_funcs[i]->F, test_apply_bndry_funcs[0]->F, Nx*Ny*NL)){
            cout << "   Function \"" << funcNames_apply_bndry[i] << "\" computes CORRECT result." << endl;
        }else{
            cout << "   Function \"" << funcNames_apply_bndry[i] << "\" computes INCORRECT result." << endl;
            exit(1);
        }
    }

    //free memory
    for(int i = 0; i<numFuncs_drift; i++){
        lbm_finalize(test_drift_funcs[i]);
    }
    for(int i = 0; i<numFuncs_set_reflective_boundaries; i++){
        lbm_finalize(test_set_reflective_boundaries_funcs[i]);
    }
    for(int i = 0; i<numFuncs_collision; i++){
        lbm_finalize(test_collision_funcs[i]);
    }
    for(int i = 0; i<numFuncs_apply_bndry; i++){
        lbm_finalize(test_apply_bndry_funcs[i]);
    }

    cout << "Program finished." << endl<<endl;
    
    return 0;
}
