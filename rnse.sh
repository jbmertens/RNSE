#! /bin/bash

rm "rnse"
gcc "RNSFluid3D-field.c" -o rnse -lm -fopenmp -lhdf5

wait

./rnse
