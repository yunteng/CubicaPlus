UTILSRC = ../../src/util/SIMPLE_PARSER.cpp \
          ../../src/util/TIMER.cpp \
          ../../src/util/TIMING_BREAKDOWN.cpp \
          ../../src/util/NNLS.cpp \
          ../../src/util/RGB_HSV.cpp

GEOMETRYSRC = ../../src/geometry/TET.cpp \
           ../../src/geometry/TRIANGLE.cpp \
           ../../src/geometry/TET_MESH.cpp \
           ../../src/geometry/TET_MESH_IO.cpp \
           ../../src/geometry/TET_MESH_INIT.cpp \
           ../../src/geometry/TET_MESH_DRAW.cpp \
           ../../src/geometry/TET_MESH_GEODESIC.cpp \
           ../../src/geometry/SUBSPACE_TET_MESH.cpp \
           ../../src/geometry/VERTEX_BVH.cpp \
           ../../src/geometry/TET_BVH.cpp \
           ../../src/dtgrid/SPARSE_SDF.cpp \
           ../../src/deformCD/aabb.cpp \
           ../../src/deformCD/aap.cpp

MATERIAL = ../../src/material/COROTATION.cpp \
           ../../src/material/FULLSPACE_COROTATION_CACHE.cpp \
           ../../src/material/SUBSPACE_COROTATION_CACHE.cpp \
           ../../src/material/PARTITIONED_SUBSPACE_COROTATION_CACHE.cpp

GLVU_LIB = /pixar/ws/users/yunteng/lib/glvu/Build/lib/Linux-gnu-opt/
GLVU_INCLUDE = /pixar/ws/users/yunteng/lib/glvu/Build/include/
GLVUFLAGS = -L${GLVU_LIB} -lglvuglvu -lglvufileutils -lglvucamera -lglvumath 
GLFLAGS = -lGL -lGLU -lglut 


CPPFLAGS = -O2 -c -fpermissive -Wno-deprecated-declarations -DUSING_LINUX -fopenmp \
            -I${GLVU_INCLUDE} \
            -I../../include/dtgrid \
            -I../../include \
            -I../../
      
LDFLAGS_COMMON = ${GLVUFLAGS} ${GLFLAGS} -fopenmp

IPOPTINCL = `PKG_CONFIG_PATH=/pixar/ws/users/yunteng/lib/ipopt/lib64/pkgconfig:/pixar/ws/users/yunteng/lib/ipopt/lib/pkgconfig:/pixar/ws/users/yunteng/lib/ipopt/share/pkgconfig: pkg-config --cflags ipopt`

IPOPTLIBS = `PKG_CONFIG_PATH=/pixar/ws/users/yunteng/lib/ipopt/lib64/pkgconfig:/pixar/ws/users/yunteng/lib/ipopt/lib/pkgconfig:/pixar/ws/users/yunteng/lib/ipopt/share/pkgconfig: pkg-config --libs ipopt`

IPOPTLINK =  -Wl,--rpath -Wl,/pixar/ws/users/yunteng/lib/ipopt/lib
