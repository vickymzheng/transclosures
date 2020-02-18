# transclosures

This repository is for reproducibility purposes, we provide all details regarding data preparation here. 
To validate our proposed approach we implemented our transitive closure and incremental component counting algorithm. 

The simulated dataset can be found here: https://drive.google.com/file/d/1itq6n2P_4MUXDMs_v6TYe3EjTUWhZ7x0/view?usp=sharing
The real world dataset can be found here: https://github.com/LomanLab/mockcommunity

Once you have downloaded the simulated data, uncompress it with the command: 
gunzip sim_reads25k.fasta.gz

Our software is located in stream_cc.cpp. In order to compile it, you must have: 
boost (to be able to use lboost_program_options)
gcc

We have tested this with boost version 1.71.0 and gcc version 4.2.1

You can compile the application with the command "make"

An example edge file from the simulated dataset can be found here: https://drive.google.com/file/d/1osPUDSJ9BBwEoegCvom9yaBM9iZ-BCfA/view?usp=sharing

You will use nano-tools-stream.cpp to stream the edge file. In order to stream the file, we use a parallel hash map. You will need this too if you would like to stream the edge file. parallel-hashmap can be found here: https://github.com/greg7mdp/parallel-hashmap
It is also included in this repository as a submodule. 

