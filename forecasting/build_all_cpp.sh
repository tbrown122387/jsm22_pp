#!/bin/bash

cd forecasting
mkdir build
cd build
cmake ..
make
touch forecasting_up_to_date.txt
