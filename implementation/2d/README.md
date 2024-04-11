
# Lattice Boltzmann Method Simulation

This project implements the Lattice Boltzmann Method (LBM) for fluid dynamics simulations in C. The primary goal of this project is to optimize a baseline python implementation of the LBM that can be used for various computational fluid dynamics simulations.

## Description

The project includes a working baseline defined in `lbm2d.c`, which accurately models fluid flow using LBM. (The baseline is a c implementation of the Python baseline we were given).

We have a working infrastructure where we can introduce optimizations to the simulation by adding functions with the same signature as the baseline to the file `lbm2d.c` and registering them for automatic execution and timing.

The infrastructure is based on the infrastructure that was used for our ASL Homeworks. 

There are two main executables provided:

- `./validate_baseline`: Verifies the consistency of the LBM baseline implementation in `lbm2d.c` with the Python-based baseline in `latticeboltzmann.py`.

- `./lbm_simulation`: Runs the baseline simulation implemented in `lbm2d.c`, executes all registered functions, tests them for correctness, and benchmarks their performance.

### Adding Optimizations

To optimize the baseline, follow these steps:

1. Implement your optimized function in `lbm2d.c` with the signature:

   ```c
   void lbm_optimized(double *F, int *cylinder, int Nx, int Ny, int n_steps);

2. Register your function in register_functions() within lbm2d.c:

    ```c
    void register_functions() {
        // Baseline runs automatically
        // Register additional functions here
        add_function(lbm_optimized, "Description of your optimization", 1);
    }



## Getting Started

### Dependencies

- A C compiler (like `gcc` or `clang`)
- CMake (minimum version 3.10)
- Python (if running the Python baseline for comparison)

### Building the Project

1. Navigate to the `2d` directory where the `CMakeLists.txt` is located:

   ```sh
   cd 2d
   ```

2. Create a new directory for the build files and navigate into it:

   ```sh
   mkdir -p build && cd build
   ```

3. Run CMake to generate the Makefiles and build the project:

   ```sh
   cmake ..
   make
   ```

This process compiles the code and generates two executables: `validate_baseline` and `lbm_simulation`.

### Running the Simulations

- To run the baseline validation, execute the following:

  ```sh
  ./validate_baseline
  ```

  This will check if the C implementation matches the results of the Python simulation.

- To run the LBM simulation with performance testing and timing, execute:

  ```sh
  ./lbm_simulation
  ```

  The program will run the baseline simulation, additional registered functions, and output performance metrics.


## Acknowledgments

- The baseline version of the algorithm is from  [this github repository](https://github.com/pmocz/latticeboltzmann-python/blob/main/latticeboltzmann.py). 
