SHELL := /bin/bash -e

BINS = ConvertTetGen RigMesh FullSim TransformDisplacements ComputeBasis PartitionedInternalForceCubature SCFCubature PartitionedHybridSim

SRCS = cubature deformCD dtgrid geometry glvu linearalgebra material util

all: 
	-for d in $(BINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -j12;cd ../..); done

cleanall: 
	-for d in $(BINS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./projects/$$d; make clean;cd ../..); done 
	-for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm -f *.o; cd ../..); done

objclean:
	-for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm -f *.o; cd ../..); done

ctags:
	ctags -R *
