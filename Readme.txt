This repository contains the code for all algorithms discussed in the paper "Mathematical models for the kidney exchange problem with reserve arcs" by Maxence Delorme, Wendy Liu, and David Manlove. 

Our algorithms are coded in C++ and use the commercial solver Gurobi for the ILP models. 
This repository containes the 16 solutions methods described in Section 7.1. 
In folder CHAIN, the 2 solution methods able to deal with non-directed donors are provided, together with the newly generated KEP instances with 10% of non-directed donors. 

The different folders correspond to the following methods in our paper:
- 1_CYCLE			| CF+NONE
- 1_CYCLE2			| CF+1RA
- 1_CYCLE3			| CF+PICORA
- 1_CYCLE3LP			| CF+PICORA+RCVF
- 3_EEF				| EEF+NONE
- 3_EEF2			| EEF+1RA
- 3_EEF3			| EEF+PICORA
- 3_EEF3LP			| EEF+PICORA+RCVF
- 4_PIEF			| PIEF+NONE
- 4_PIEF2			| PIEF+1RA
- 4_PIEF3			| PIEF+PICORA
- 4_PIEF3LP			| PIEF+PICORA+RCVF
- 5_HCF				| HCF+NONE
- 5_HCF2			| HCF+1RA
- 5_HCF3			| HCF+PICORA
- 5_HCF3LP			| HCF+PICORA+RCVF
- CHAINS/1_CYCLE3LP_CHAIN1	| CF+PICORA+DPICEF+RCVF
- CHAINS/1_CYCLE3LP_CHAIN3	| CF+PICORA+IPICEFTH4+RCVF

Each folder contains the same substructure. For example, 1_CYCLE contains the following files:
- Allocation.cpp		| Contains a number of secondary functions (this file is usually the same for each subfolder)
- Allocation.h			| The header file corresponding to Allocation.cpp (this file is usually the same for each subfolder)
- main.cpp			| The front-end code for using the method  
- main.h			| The header file corresponding to main.cpp 
- makefile			| Used for compiling under linux (it needs to be updated by the user)
- time.cpp			| A generic file used to measure the computation time 
- time.h			| The header file corresponding to time.cpp 

********

Once compiled, the following command can be used to run the algorithm:
	./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" "K" "B" for the 16 first approaches
	./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" "K" "L" "B" for the 2 last approaches
where
- PROGRAM is the name of the compiled software 
- ./PATH_INSTANCE is the relative path of the folder where the instance to solve is located
- NAME_INSTANCE is the name of the instance to solve
- ./PATH_AND_NAME_OUTPUT_GENERAL is the name of the file (together with its relative path) where performance metrics (such as the optimality status, the CPU time required, or the number of variables) are stored after solving an instance
- K is the maximum cycle size
- L is the maximum chain length (only for the models allowing chains)
- B is the maximum number of reserve arcs allowed 

********

Moreover, "_INPUT.rar" contains a txt-file for the nexly generated test instances. 
See "https://wpettersson.github.io/kidney-webapp/#/" for a detailed explanation of the generated instance files.

