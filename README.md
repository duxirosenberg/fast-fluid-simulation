
# ETHZ Advanced Systems Lab 24 - Team 54
# Fast Fluid Simulation: Optimizing the Lattice Boltzmann Method

## Overview

This repository contains the code and documentation for our project on optimizing the Lattice Boltzmann Method (LBM) for fluid dynamics simulations. Our primary goal was to enhance the performance of LBM by implementing various optimizations, resulting in a significant speedup.

For more details on our results, refer to the full project report included in this repository [TODO: Link].

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Optimization Techniques](#optimization-techniques)
- [Experimental Results](#experimental-results)
- [Contributors](#contributors)
- [References](#references)
- [Notes on Infrastructure](#notes-on-infrastructure)

## Introduction

Fluid dynamics simulations are vital in engineering and environmental sciences. Traditional methods like solving the Navier-Stokes equations can be computationally intensive. The Lattice Boltzmann Method (LBM) provides a more efficient alternative by solving the Boltzmann transport equation on a discrete lattice grid. This project presents optimized implementations of LBM, focusing on three boundary conditions: Lees-Edwards, Couette, and periodic, achieving a speedup of up to 90x.

We used a C++ implementation from Callum Marshall as a starting point, which was heavily optimized for performance. Our work built upon this foundation, adapting it into C and applying further optimizations.

## Features

- **Optimized LBM Implementation:** Significantly improved performance for fluid dynamics simulations.
- **Boundary Conditions:** Support for Lees-Edwards, Couette, and periodic boundary conditions.
- **Extensive Testing Infrastructure:** Comprehensive testing and timing setup to validate performance gains.
- **Optimization Techniques:** Incorporates strength reduction, loop ordering, vectorizations, and more.

## Installation

### Prerequisites

- GCC with support for AVX2 instructions
- Python (for plotting and analysis scripts)

### Steps

1. Clone the repository:

2. Create a build directory, and compile the Code:

    ```
    mkdir build && cd build
    cmake .. && make
    ```

3. Install Python dependencies for plotting (optional):
    ```sh
    pip install matplotlib pandas
    ```


## Usage

### Running the Simulation

This repository provides two main executables: `timing` and `cmdline`.

#### Timing

The `timing` executable is used for testing and timing various parts of the simulation. It can be executed with the following parameters:
```sh
./timing $DIRECTION $NX $NY $NZ $STREAM $MOMENTUM $COLLISION $LBM $RESET_DATAFILE
```
- `$DIRECTION`: Direction of the simulation.
- `$NX`, `$NY`, `$NZ`: Grid sizes in the x, y, and z dimensions, respectively.
- `$STREAM`, `$MOMENTUM`, `$COLLISION`,`$LBM`: Flags to specify which parts of the simulation should be tested and timed.
- `$RESET_DATAFILE`: Flag to reset the data file.

To run the timing executable for various grid sizes and produce all necessary experimental results, you can use the `run.sh` script located in the `visualization` directory. Note that executing this script may take a long time.

#### Cmdline

The `cmdline` executable is the main entry point for running the full LBM simulation based on the configurations specified in the `options.json` file. Here’s a quick overview of what it does:

- Checks for the existence of `options.json`.
- Parses simulation parameters such as grid size, velocity set, boundary conditions, number of steps, and more.
- Initializes arrays and structures for the simulation.
- Optionally generates initial conditions.
- Loads initial conditions if `ic.csv` exists; otherwise, uses default values.
- Runs the LBM simulation for the specified number of steps, saving data periodically based on the configuration.

To run the `cmdline` executable:
```sh
./cmdline
```

Ensure `options.json` is correctly configured with the desired simulation parameters. The program will prompt for certain inputs and output results to the specified directories.

### Generating Plots

We used Jupyter notebooks, which can be found in the `visualization` directory to generate roofline, performance and runtime plots.



## Optimization Techniques

1. **Strength Reduction and Precomputing Loop Invariants:** Simplifying arithmetic operations to reduce computational load.
2. **Cache Optimizations:** Reordering loops to improve data locality and reduce cache misses.
3. **Memory Layout Adaptations and `memcpy`:** Using efficient memory operations to optimize data movement.
4. **Manual Vectorization with AVX2:** Leveraging advanced vector instructions to accelerate floating-point operations.

## Experimental Results

Our optimizations resulted in substantial performance improvements, especially for larger grid sizes. The highest speedup achieved was 90x for the Couette boundary condition at a grid size of 64³. Detailed results and roofline plots are available in the report and in the **`visualization`** directory.

## Contributors

- **Johannes Gasser:** Core implementation, timing infrastructure, periodic and Lees-Edwards boundary condition optimization.
- **Severin Klapproth:** Struct-based implementation, executables, collision function optimization, further optimization experiments
- **Douglas Orsini-Rosenberg:** 2D to C conversion, timing infrastrucutre, Couette streaming optimization.
- **Rachel Schuchert:** Boundary condition comparison, periodic and momentum function optimization, and all results visualization/plotting.

## References

1. Aidun, C. K., & Clausen, J. R. (2010). Lattice Boltzmann Method for Complex Flow Simulations. *Annual Review of Fluid Mechanics*, 42(1), 439–472.
2. Krüger, T., Kusumaatmaja, H., Kuzmin, A., Shardt, O., Silva, G., & Viggen, E. M. (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer.
3. Mocz, P. (2020). Create your own lattice Boltzmann simulation (with Python). Princeton University. [Link](https://medium.com/swlh/create-your-own-lattice-boltzmann-simulation-with-python-8759e8b53b1c).
4. Marshall, C. (2019). LBM. GitHub. [Link](https://github.com/callummarshall9/LBM).
5. Smith, Z. (2012). Bandwidth: a memory bandwidth benchmark. [Link](https://zsmith.co/bandwidth.php).

For more details, refer to the full project report included in this repository.

## Notes on Infrastructure

### Directory Structure

- **`src` directory:** Contains the actual C code for the simulation.
- **`include` directory:** Contains all the headers.
- **`scripts` directory:** Contains the infrastructure and scripts for performing timing and testing.
- **`visualization` directory:** Contains Python scripts, notebooks, and tools for plotting the experimental results.
- **`options.json`:** Configuration file to change some of the simulation setup.

### Adding New Optimizations

To optimize a part of the simulation, follow these steps:

1. Navigate to the source file you want to optimize, for example, `src/3d/collision.c`.
2. Add the new function with your optimizations.
3. Update the corresponding header file in `include/3d/collision.h` and register the new function in `register_collision_functions()`.

The new function will be automatically tested and timed when running the scripts. This infrastructure is based on the setup used in the Code Expert homeworks from the course.

With this setup, we were able to efficiently develop, test, and benchmark various optimizations to the LBM simulation.
