cmake_minimum_required(VERSION 3.21)
project(jsmpp_v2)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -pthread)

include_directories(include)
include_directories(/usr/include/eigen3) # installed eigen3 here
include_directories(/usr/local/include)  # installed ssme, pf and catch2 here

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Catch2 3)
find_package(Threads REQUIRED)
find_package(ssme CONFIG REQUIRED)
find_package(pf CONFIG REQUIRED)

# list all source files explicitly
# https://stackoverflow.com/questions/1027247/is-it-better-to-specify-source-files-with-glob-or-each-file-individually-in-cmak/18538444#18538444
add_executable(jsmpp_v2 main.cpp)
