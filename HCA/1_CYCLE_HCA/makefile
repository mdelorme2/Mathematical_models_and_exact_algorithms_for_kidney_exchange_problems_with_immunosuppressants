GUROBI_DIR = /home/mdelorme/data/volume_2/gurobi1102/linux64
CLIB     = -L$(GUROBI_DIR)/lib/ -lgurobi110
INC      = $(GUROBI_DIR)/include/
CPPLIB   = -L$(GUROBI_DIR)/lib/ -lgurobi_c++ $(CLIB)
GUROBI_OPTS = -I$(INC) $(CPPLIB) -lpthread -lm -m64 

CC = g++ 
CFLAGS = -O2 -Wall -ansi -pedantic -DIL_STD -std=c++11
DEBUG = -pg -g -Wall -ansi -pedantic -DIL_STD 
OBJECTS = main.o Allocation.o time.o 
 
exec : $(OBJECTS) 
	$(CC) $(CFLAGS) -o CYCLE $(OBJECTS) $(GUROBI_OPTS)
 
.cpp.o : 
	$(CC) $(CFLAGS) $(GUROBI_OPTS) -c $< -o $@ 
 
clean : 
	rm -f $(OBJECTS) CYCLE 