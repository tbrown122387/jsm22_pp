# Posterior Predictive Distribution Project

This project is primarily organized into a bunch of `R` scripts that clean data, analyze data, and then plot output. If you'd like to see how these files relate to one another, there is a nice visualization in `file_dep_vis/file_dep_vis.svg`. 

These dependencies are also encoded in the `Makefile`, which is a utility that you should run the code from. There will be some example commands below that will allow you to execute high-level tasks. 

Before any of the analysis code can be run, the C++ code must be compiled. Sometimes that can be tricky, particularly for users who have more experience using interpreted languages such as `R`.


## C++ Code

### Third-Party Software Requirements

To build the c++ code, you need:

 1. cmake to build c++ code
 2. Eigen3 for linear algebra functionality
 3. Boost for other math functionality
 4. [`Catch2`](https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#installing-catch2-from-git-repository) for the unit tests
 5. [`pf`](https://github.com/tbrown122387/pf) a particle filtering library
 6. [`ssme`](https://github.com/tbrown122387/ssme) a more experimental particle filtering/MCMC library

### How to Build C++ Code

The cpp code is in `cpp/`. After installing all the above software, you can build the project with CMake. Most IDEs handle CMake projects, and you can also do it from the command line. 

Note: depending on where you installed the above software, you might have to chane a few lines in `cpp/CMakeLists.txt`. For example, say you downloaded `eigen` and didn't install it to `/usr/include/eigen3`. Then the line `include_directories(/usr/include/eigen3)` in `cpp/CMakeLists.txt` is wrong for you, and you need to change the path in there. 

#### Option 1: CLion IDE

The `cpp/` folder is a CMake project that can be edited with CLion. To open up the `cpp/` code as a CLion project,
[click here.](https://www.jetbrains.com/help/clion/creating-new-project-from-scratch.html) From there, just hit the build button.

#### Option 2: Command Line Build

The (Debug) build command I use to build everything on Linux/Mac is

	/path/to/cmake/bin/cmake --build /path/for/software/jsm22_pp/cpp/cmake-build-debug --target jsmpp_v2 -- -j 6

There is also a similar command to build the Release configuration.

### How to Run

#### Option 1: CLion IDE

If you build with CLion, you can run with CLion. This program requirements command line arguments, though. You need to specify the "run mode", and the path to the return observations. For help choosing the run mode, look at the comments in `cpp/main.cpp`. For help instructing CLion to use arguments when it runs the program, [click here](https://www.jetbrains.com/help/clion/run-debug-configuration.html#envvars-progargs).

#### Option 2: Command line interface

You can also run the program from your terminal. Here's an example of how to run the "Release" build. But remember, file paths may change, and you may change the "run mode."

	./cpp/cmake-build-release/jsmpp_v2 1 "data/SPY_returns.csv"

You should run this command from the root directory because there are a few file paths hardcoded in the files :( R scripts that call c++ code follow this approach--they construct a command similar to the one above, and then pass it to `system()`.

# Running Everything

Once all the c++ code is built, you can start running everything my simply typing `make` in the command line. NB: you can also run R scripts on a case-by-case basis--to do that, follow the examples in `Makefile` (e.g. `Rscript yadayadayada.R`). 

