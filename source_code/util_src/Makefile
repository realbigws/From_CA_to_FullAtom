

# make the subdir first
CXX = g++
FLAG = -O3
SHELL :=/bin/bash

all:
	#-> PDB operation
	$(CXX) PDB_To_SEQ_miss.cpp $(FLAG) -o PDB_To_SEQ_miss
	$(CXX) Protein_DynaProg.cpp $(FLAG) -o Protein_DynaProg
	$(CXX) PDB_Add_Chain.cpp $(FLAG) -o PDB_Add_Chain
	#-> move to previous folder
	mv PDB_To_SEQ_miss Protein_DynaProg PDB_Add_Chain ../

# clean all the temp files
clean:
