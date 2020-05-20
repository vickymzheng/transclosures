# transclosures

This repository is for reproducibility purposes, we provide all details regarding data preparation here. 
To validate our proposed approach we implemented our transitive closure and incremental component counting algorithm. 

The simulated dataset can be found here: https://drive.google.com/file/d/1itq6n2P_4MUXDMs_v6TYe3EjTUWhZ7x0/view?usp=sharing
The real world dataset can be found here: https://github.com/LomanLab/mockcommunity

Once you have downloaded the simulated data, you can uncompress it with the command: 
`gunzip sim_reads25k.fasta.gz`

An example edge file from the simulated dataset can be found here: https://drive.google.com/file/d/1osPUDSJ9BBwEoegCvom9yaBM9iZ-BCfA/view?usp=sharing

Edge files must follow the format:

```
read_id1	read_id2	readlen1	readlen2	overlaplen	read1_orient	read2_orient
```


Our software is located in stream_cc.cpp. In order to compile it, you must have: 
`boost` (to be able to use `lboost_program_options`) and
`gcc`

We have tested this with `boost` version 1.71.0 and `gcc`version 4.2.1

You can compile the application with the command `make`

You will use nano-tools-stream.cpp to stream the edge file. In order to stream the file, we use a parallel hash map. You will need this too if you would like to stream the edge file. parallel-hashmap can be found here: https://github.com/greg7mdp/parallel-hashmap
It is also included in this repository as a submodule. 

You can stream an edge file with the command:

```./nano-tools-stream [edge file] [random seed]```

With the edge file provided, you can stream it with the command:

 ```./nano-tools-stream sim-reads25k_ov1000.edges 1```

If you want to feed this stream into our program for processing, you can do this using the pipe like so:

 ```./nano-tools-stream sim-reads25k_ov1000.edges 1 | ./stream_cc```

You can get instructions on how to use our program with the command:

`./stream_cc -h`

stream_cc produces 6 output files: 
* default_out.component_convergence - Each line entry is the number of components seen at that time. So if you see the number 8 on line 199990, this means there are 8 components after processing the 199990th node.
* default_out.component_sizes - The first column of this file gives component ids while the second column gives the size. This lists all components found at the end of processing
* default_out.irreducible_edges - This file gives our final irreducible graph at the end of processing. The first column is source nodes and the second column is target nodes.
* default_out.storage_rate - Each line entry is the number of nodes stored at that time. So if you see the number 3010 on line 199990, this means there are 3010 nodes stored after processing the 199990th node.
* default_out.trans_mapping - This gives the component mapping of each node (from the paper, this is UF). The first column is a node id and the second column is the component it is mapped to.
* default_out.trans_stats - This gives a summary of our application's findings.

Note that on the example data set provided, our application took about 2 minutes to finish running on 200,00 reads. 


