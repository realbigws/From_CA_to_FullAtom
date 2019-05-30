

# make the subdir first
CXX = g++
FLAG = -O3
SHELL :=/bin/bash

all:
	#-> build3Dmodel
	$(CXX) build3Dmodel.cpp $(FLAG) -o build3Dmodel
	$(CXX) build3Dmodel_long_distance.cpp $(FLAG) -o build3Dmodel_long_distance
	mkdir -p ../buildModel
	mv build3Dmodel build3Dmodel_long_distance ../buildModel
	#-> util
	cd util_src&&make
	mkdir -p ../util
	mv PDB_To_SEQ_miss Protein_DynaProg PDB_Add_Chain ../util

# clean all the temp files
clean:
