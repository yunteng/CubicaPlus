UTILSRC = ../../src/util/SIMPLE_PARSER.cpp \
					../../src/util/TIMER.cpp \
					../../src/util/TIMING_BREAKDOWN.cpp \
					../../src/util/NNLS.cpp \
					../../src/util/RGB_HSV.cpp \
					../../src/glvu/glvu.cpp \
					../../src/linearalgebra/BLOCK_COO_MATRIX.cpp \
					../../src/linearalgebra/BLOCK_SPARSE_MATRIX.cpp

GEOMETRYSRC = ../../src/geometry/TET.cpp \
					 ../../src/geometry/TRIANGLE.cpp \
					 ../../src/geometry/TET_MESH.cpp \
					 ../../src/geometry/TET_MESH_IO.cpp \
					 ../../src/geometry/TET_MESH_INIT.cpp \
					 ../../src/geometry/TET_MESH_DRAW.cpp \
					 ../../src/geometry/SUBSPACE_TET_MESH.cpp \
					 ../../src/geometry/VERTEX_BVH.cpp \
					 ../../src/geometry/TET_BVH.cpp \
					 ../../src/dtgrid/SPARSE_SDF.cpp \
					 ../../src/deformCD/aabb.cpp \
					 ../../src/deformCD/aap.cpp \
					 ../../src/geometry/SURFACE.cpp \
					 ../../src/geometry/SPHERE.cpp \
					 ../../src/geometry/CYLINDER.cpp \
					 ../../src/geometry/PLANE.cpp

MATERIAL = ../../src/material/COROTATION.cpp \
					 ../../src/material/FULLSPACE_COROTATION_CACHE.cpp \
					 ../../src/material/PARTITIONED_FULLSPACE_COROTATION_CACHE.cpp \
					 ../../src/material/SUBSPACE_COROTATION_CACHE.cpp \
					 ../../src/material/PARTITIONED_SUBSPACE_COROTATION_CACHE.cpp

CUBATURE = ../../src/cubature/SCF_CUBATURE_LOADER.cpp

GLFLAGS = -L/usr/lib/ -lGL -lGLU -lglut

CC = g++

CPPFLAGS = -m64 -O2 -fpermissive -Wno-deprecated-declarations -fopenmp -c \
						-I../../ \
						-I../../include \
						-I../../include/dtgrid \
						-I../../include/glvu
 
LDFLAGS_COMMON = ${GLFLAGS} -fopenmp
