#!/bin/bash

NUM_MCMC_ITERS=10
NUM_PFILTERS=1

# backup old samples and messages files
echo "backing up old samples and messages files"
find $PWD/data -type f -name 'samps_*' -not -name "*.bck" -print0 | xargs --null -I{} mv {} {}.bck
find $PWD/data -type f -name 'messages_*' -print0 | xargs --null -I{} mv {} {}.bck

# rebuild estimation program
echo "rebuilding estimation program"
cd estimation
mkdir build
cd build
cmake ..
make

# run estimation program
echo "sampling from the posterior"
cd ..
./build/svol_leverage_estimation ../data/SPY_returns.csv samps messages $NUM_MCMC_ITERS $NUM_PFILTERS


# copy recent parameter samples file and move it to param_samples.csv
cp messages_* ../data
cp samps_* ../data
find . -type f -name "samps_*" -not -name "*.bck" -exec cp "{}" ../data/param_samples.csv \;


