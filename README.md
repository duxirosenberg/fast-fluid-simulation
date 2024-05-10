# Building
Create a build directory in this directory & build with `cmake`:
```
mkdir build
cd build
cmake .. && make
```

# Running
Add a valid `options.json` to the build directory. There should exist an example in the main directory of the 3d implementation.
In the build directory run the `main_cmdline` executable. 

# Timing
In the build directory run the `main_timing` executable

## Adding more functions
Create a C file in `clib/source` containing the function either using the struct or the array approach.
In the `LBMFunctions.h` header, add your function to the header and register your function in the 
`register_functions` method. Build the project again and run the timing executable again to see your functions results.
