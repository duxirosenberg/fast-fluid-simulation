
# Lattice Boltzmann Method Simulation

This project implements the Lattice Boltzmann Method (LBM) for fluid dynamics simulations in C. The primary goal of this project is to optimize a baseline C++ implementation of the LBM that can be used for various computational fluid dynamics simulations.

## Description

The project includes a working baseline defined in `LBMarrays.c` (TODO CURRENTLY THERE IS ALSO LBMbaseline.c) and `LBMstructs.c`, which contains the straighforward C baselines one using different arrays and one using a sutrct to store all the data. (The baseline is a c implementation of the C++ baseline we were given).

We have a working infrastructure where we can introduce optimizations to the simulation by adding functions with the same signature as the baseline.

The infrastructure is based on the infrastructure that was used for our ASL Homeworks. 

There are two main executables provided:

- `./cmdline`: executable that runs the current implementation using predefined options.

- `./timing`: executable that runs baseline and all registered functions for various problem sizes. The results can be saved to csv. Additionally, the correctness of each registered function is tested against the baseline.

There is also the `plots.ipynb` notebook to generate performance and roofline plots.

### Adding Optimizations

To optimize the baseline:

Create a C file in `clib/source` containing the function either using the struct or the array approach.
In the LBMFunctions.h header, add your function to the header and register your function in the register_functions method. Build the project again and run the timing executable again to see your functions results.



## Getting Started

### Dependencies

- A C compiler (like `gcc` or `clang`)
- CMake (minimum version 3.10)


### Building the Project

1. Navigate to the `3d` directory where the `CMakeLists.txt` is located:

   ```sh
   cd 3d
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



### Running the Simulations
Add a valid options.json to the build directory. There should exist an example in the main directory of the 3d implementation.
In the build directory run the main_cmdline executable.

In the build directory run the main_timing executable


## Acknowledgments

- The baseline version of the algorithm is from  [this github repository](https://github.com/callummarshall9/LBM/tree/master). 
