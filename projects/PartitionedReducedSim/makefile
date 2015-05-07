include ../common.inc

ODE = ../../bin/PartitionedReducedSim

SOURCEC    = $(UTILSRC) \
						 $(GEOMETRYSRC) \
						 $(MATERIAL) \
						 $(CUBATURE)

OBJSC = $(SOURCEC:.cpp=.o)


.cpp.o:
	$(CC) $(CPPFLAGS) $< -o $@
	
regression: ALL
ALL: $(SOURCES) ode

ode : ode.o $(OBJSC)
	${CC} ode.o $(OBJSC) ${LDFLAGS_COMMON} -o ${ODE}

clean:
	rm -rf *.o ${ODE}

cleanall:
	rm -rf *.o $(OBJSC) $(ODE)

