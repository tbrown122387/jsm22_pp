# Posterior Predictive Distribution Project

This project is organized into a few parts. The most complex part is all of the C++ code. This handles all of the computationally intensive tasks. The rest is a bunch of scripts that either call this code, edit config files for it, or work with output from it. 

## C++ Code

### Third-Party Software Requirements

You need:

 1. cmake to build c++ code
 2. Eigen3 for linear algebra stuff
 3. Boost for math stuff
 4. [`Catch2`](https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#installing-catch2-from-git-repository) for the unit tests
 5. [`pf`](https://github.com/tbrown122387/pf) a particle filtering library
 6. [`ssme`](https://github.com/tbrown122387/ssme) a more experimental particle filtering/MCMC library

### How to Build C++ Code

The cpp code is in `cpp/`. After installing all the above software, you can build the project with CMake. Most IDEs handle CMake projects, and you can also do it from the command line.

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

# Everything Else

TODO

