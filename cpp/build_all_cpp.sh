#!/bin/bash

cd ~/jsm22_pp/cpp
mkdir build
cd build
cmake ..
make
touch cpp_up_to_date.txt
