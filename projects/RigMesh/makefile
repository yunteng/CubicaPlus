include ../common.inc

ODE = ../../bin/RigMesh

SOURCEC    = $(UTILSRC) \
             $(GEOMETRYSRC) \
             $(MATERIAL)

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

