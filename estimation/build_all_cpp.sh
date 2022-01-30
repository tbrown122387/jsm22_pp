#!/bin/bash

cd estimation
mkdir build
cd build
cmake ..
make
touch estimation_up_to_date.txt
