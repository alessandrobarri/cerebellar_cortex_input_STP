
# Synaptic basis of a sub-second representation of time

This repository contains the c++ code used for carying out the simulations of the publication: "Synaptic basis of a sub-second representation of time". The code 1) generates instances of cerebellar-cortex networks with the desired statistical properties, as described in the publication (i.e. mossy fiber rates, release probablities etc.); 2) carries out the numerical integration of the dynamics for the specified network(s); and 3) performs learning of the granule cell to purkinje cell weights for the eye-lid conditioning or Bayesian interval estimation task.

## System and Software requirements

The code is written in c++ and requires the armadillo library (http://arma.sourceforge.net/) and the GNU scientific library (gsl, https://www.gnu.org/software/gsl/).

Compilation and excecution of the code has been tested on Ubuntu 20.04. with gcc/g++ version 9.4.0., gsl version 2.6 and armadillo version 11.1.

No non-standard hardware is required to run the program.

## Installation guide

To install the program, copy all file in a custom folder. The code comes with two folders that need to be in the main folder in order for the program to work correctly: 1) a folder named 'headers' which contains essential sub-routines; 2) a folder named 'output' into which the program will write the simulated data.

There are two main program files: 'simCC_learn-MULTI_2pools_full.cpp' and 'simCC_learn-MULTI_2pools.cpp', which simulate the full version (Figures 1,2 and 6) and the reduced version (Figures 3,4,5 and 6) of the cerebellar cortex model, respectively. Under linux, the code can be compiled using:

g++ -O3 -std=c++11 -o simCC_learn-MULTI_2pools.out simCC_learn-MULTI_2pools.cpp -lgsl -lgslcblas -lm -larmadillo

g++ -O3 -std=c++11 -o simCC_learn-MULTI_2pools_full.out simCC_learn-MULTI_2pools_full.cpp -lgsl -lgslcblas -lm -larmadillo

Here, -lgsl -lgslcblas -lm and -larmadillo are necessary to link the gsl and armadillo libraries, respactively. Note that the minimum required C++ standard for running armadillo is c++11.

Copying the files and compiling the program should not take more than 5 minutes overall.

## Demo


In the 'sample' mode, the programs should run between 40s to 90s (for 4000 learning iterations).

## example input/output (in example folder)

## run on ubuntu, approximate time for execution

## license?

