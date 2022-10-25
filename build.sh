#!/bin/sh -e

#CPPFLAGS="-std=c++20 -ggdb3 -march=native -Wall -Wextra -Wpedantic"
CPPFLAGS="-std=c++20 -O3 -march=native -fopenmp"
EIGEN_INC="./eigen"

c++ $CPPFLAGS -I"$EIGEN_INC" main.cpp -o main