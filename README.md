# Lattice Boltzmann Method Simulation

This project implements the Lattice Boltzmann Method (LBM) for fluid dynamics simulations in C. The primary goal of this project is to optimize a baseline C++ implementation of the LBM that can be used for various computational fluid dynamics simulations.

## Description
This project implements a c baseline of the LBM method. The LBM simulation performs 3 consecutive steps on a particle distribution:


-

The project includes a working baseline which contains the straighforward C baselines one using different arrays and one using a sutrct to store all the data. (The baseline is a c implementation of the C++ baseline we were given).

We have a working infrastructure where we can introduce optimizations to the simulation by adding functions with the same signature as the baseline.

The infrastructure is based on the infrastructure that was used for our ASL Homeworks. 

There are two main executables provided:

- `./cmdline`: executable that runs the current implementation using predefined options.

- `./timing`: executable that runs baseline and all registered functions for various problem sizes. The results can be saved to csv. Additionally, the correctness of each registered function is tested against the baseline.

There is also the `plots.ipynb` notebook to generate performance and roofline plots.


# How to work with this code

## 1. Building
To create a build directory, run the following commands:

```
mkdir build && cd build
cmake .. && make
```

This creates several executables:

- `./cmdline`:
- `./timing`:

## 2. Running
Run either of the 2 executables
- `./cmdline`:
- `./timing`:


# Timing
In the build directory run the `./timing` executable

# Optimizations

## Adding more functions
If you want to optimize a step of the simulation, go to the corresponding c file in 
src/3d directory and add a new function. 

Head to the correspondign header and register your function using the `register_functions` method. Build the project again and run the timing executable again to see your functions results. 


## Acknowledgments

The baseline version of the algorithm is an adaptation of the code from  [this github repository](https://github.com/callummarshall9/LBM/tree/master). 
