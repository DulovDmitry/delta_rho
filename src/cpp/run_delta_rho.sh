#!/bin/bash

g++ -O3 -fopenmp main.cpp Grid.cpp Molecule.cpp Atom.cpp BasisFunction.cpp Orbital.cpp -lopenblas -o main -Wall -std=c++20
./main
