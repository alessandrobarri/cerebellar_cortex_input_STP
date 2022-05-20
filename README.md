
# Synaptic basis of a sub-second representation of time

This repository contains the c++ code used for carying out the simulations of the publication: "Synaptic basis of a sub-second representation of time". The code 1) generates instances of cerebellar-cortex networks with the desired statistical properties, as described in the publication (i.e. mossy fiber rates, release probablities etc.); 2) carries out the numerical integration of the dynamics for the specified network(s); and 3) performs learning of the granule cell to purkinje cell weights for the eye-lid conditioning or Bayesian interval estimation task.

## System and Software requirements

The code is written in c++ and requires the armadillo library (http://arma.sourceforge.net/) and the GNU scientific library (gsl, https://www.gnu.org/software/gsl/).

Compilation and excecution of the code has been tested on Ubuntu 20.04. with gcc/g++ version 9.4.0., gsl version 2.6 and armadillo version 11.1.

Two files are included for importing the program's output into matlab (version R2019b).

No non-standard hardware is required to run the program.

## Installation guide

To install the program, copy all files in a custom folder. The code comes with two folders that need to be in the main folder in order for the program to work correctly: 1) a folder named 'headers' which contains essential sub-routines; 2) a folder named 'output' into which the program will write the simulated data.

There are two main program files: `simCC_full.cpp` and `simCC.cpp`, which simulate the full version (Figures 1,2 and 6) and the reduced version (Figures 3,4,5 and 6) of the cerebellar cortex model, respectively. Under linux, the code can be compiled using:

` g++ -O3 -std=c++11 -o simCC.out simCC.cpp -lgsl -lgslcblas -lm -larmadillo`

` g++ -O3 -std=c++11 -o simCC_full.out simCC_full.cpp -lgsl -lgslcblas -lm -larmadillo`

Here, -lgsl -lgslcblas -lm and -larmadillo are necessary to link the gsl and armadillo libraries, respectively. Note that the minimum required C++ standard for running armadillo is c++11.

Copying the files and compiling the program should not take more than 5 minutes overall.

## Running the program and producing test output

Under linux, the code can be executed via the command line. It is necessary to respect the following syntax:

` ./simCC.out \[OUTPUT FOLDER] \[MODE] \[NP] \[NR] \[NS] \[RANDOM NUMBER SEED]`

` ./simCC_full.out \[OUTPUT FOLDER] \[MODE] \[NP] \[NR] \[NS] \[RANDOM NUMBER SEED]`

The command line arguments are:

* OUTPUT FOLDER: The sub-folder into which the output is written. Do put a slash (\/) at the end of the folder name.
* MODE: The operation mode of the program. There are three options for the full model (*cI*, *bayesian* and *sample*) and six for the reduced model (*MF_U_corr*, *Ntrials* and *scan*, in addition to the other three).
	- *sample*: runs NR instances of eye-lid conditioning learning
	- *bayesian*: runs NR instances of Bayesian interval learning
	- *cI*: runs NR instances of eye-lid conditioning learning while varying the parameter JI over NP values
	- *Ntrials*: runs NR instances of eye-lid conditioning learning while varying the number of learning iterations
	- *MF_U_corr*: runs NR instances of eye-lid conditioning learning while varying the rank-correlation between MF rates and release probabilies over NP values
	- *scan*: carries out a scan over two MF rate parameters (two out of &mu;<sub>D</sub>, &mu;<sub>S</sub>, &sigma;<sub>D</sub> ,&sigma;<sub>S</sub>) over an NPxNP grid by running NR instances of eye-lid conditioning learning
* NP: The number of parameters over which the program iterates when it is in the *cI* or *MF_U_corr* mode. In the *scan* mode, it will scan over a square grid with side length NR.
* NR: The number of different random network runs per parameter set.
* NS: The number of iterations.
* RANDOM NUMBER SEED: Enter an integer random number seed.

To test that the program works correctly, execute the following for the reduced and full model, respectively:

` ./simCC.out test/ sample 1 1 4000 -3 `

` ./simCC_full.out test_full/ sample 1 1 4000 -3 `

Upon execution, the shell prompts a line stating which operation mode has been chosen, as well as a few network parameters. Either program first creates and subsequently writes into the sub-folders `output/test/` or `output/test_full/`. If everything works correctly, the resulting simulated data should be of the same format as in the example output folders provided with the code.

In the *sample* mode, the programs should run between 40s to 90s (for NS=4000 and NR=1).

## Example output



To facilitate the handling of the program's output we included the files `read_data_simCC.m` and `read_data_simCC_bayes.m` to the repository. These scripts can be used to import simulated data and simulation parameters into matlab and organise them in an easily readable structure.

