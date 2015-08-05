============================================
Introduction

This is Cubica++, a toolkit for efficient finite element simulations of
articulated deformable bodies containing both geometric and material non-
linearities. In addition to providing the core features of Cubica
(http://www.mat.ucsb.edu/~kim/cubica/), it implements the methods described in
the papers

Simulating articulated subspace self-contact.
  Yun Teng, Miguel A. Otaduy, and Theodore Kim. SIGGRAPH 2014
Subspace Condensation: Full Space Adaptivity for Subspace Deformations
  Yun Teng, Mark Meyer, Tony DeRose, and Theodore Kim. SIGGRAPH 2015

Various versions of this code were used to generate the examples in those
papers.

============================================
Third-party libraries

To minimize the effort to compile the code, aLl necessary third-party libraries are included in this release. 

Cubica++ uses Eigen (http://eigen.tuxfamily.org/) for the linear algebra 
routines. I did a slight in Eigen/src/SparseCore/SparseUtil.h: the following 
three lines were added in the Triplet class:

Index& row() { return m_row; }
Index& col() { return m_col; }
Scalar& value() { return m_value; }

Other third party libraries include:
  dtgrid(https://code.google.com/p/dt-grid/), for generating as well as 
    looking up narrow-banded signed distance field (SDF).
  deformCD(http://gamma.cs.unc.edu/DEFORMCD/), for BVH triangle-triangle
    collision.
  glvu(http://www.cs.unc.edu/~walk/software/glvu/), for OpenGL navigation, etc.
    All calls to GLVU are restricted to include/util/VIEWER.h and
    include/util/VIEWER.inl so you can easily swap to your own favoriate 
    OpenGL navigator.

============================================
Compiling the code

Under projects: replace common.inc with one of the following: 
common.osx.noopenmp.inc, common.osx.openmp.inc, common.linux.inc. 
The osx.openmp one uses OpenMP/Clang (https://clang-omp.github.io/) which I 
recall took some effort to install. If you don't want to bother with getting
OpenMP to work, simply set USING_OPENMP and USING_SUBSPACE_OPENMP to 0 in
include/SETTINGS.h.

Once all the necessary changes are made in common.inc, call make from the
cubica++ directory. It will recurse down into the projects directories
and build all the necessary binaries. The binaries will be deposited in the bin directory.

The pipeline.sh script executes the entire workflow on the hand example, including full space simulation training, basis generation, cubature training and subspace condensation simulation.

============================================
Usage:

All binaries uses our customized config file. Please see the config direcotry 
for a few examples. Specifically, chand.fullsim.%d.cfg are for full space 
simulation, chand.training.cfg trains the subspace and chand.run.cfg is for simulation with subspace condensation.

From cubica++ directory, execute the following commands:

./bin/ConvertTetGen ./cfg/chand.fullsim.0.cfg
This call converts TetGen mesh to cubica++ format. You need to have an [tet mesh 
name].ele and [tet mesh name].node file in the "output" directory.

./bin/RigMesh ./cfg/chand.fullsim.0.cfg
This call first constrains the tets that are penetrated by the skeleton and 
writes out a new constrained mesh. Then it reloads the contrained mesh and 
computes the volume diffusion skinning weights.

./fullsim.sh
This calls ./bin/FullSim on 8 chand.fullsim config files. Each config file
points to an individual joint exerciese sequence, thumb, index, etc.. Self-
collision detection are enabled. The mesh deformation and self-collision
vertices as well as collision responses are saved to data path. This stage
takes about 1 hour to finish. You can run config file 0-6 in parallel to
accelerate the process. chand.fullsim.7.config requires starting from frame 
148 which is computed in chand.fullsim.6.config.

./bin/TransformDisplacements ./cfg/chand.training.cfg
This transforms the world-space displacement vectors to before-skinning
displacements. They are used to compute the skinning-correction basis. You
can skip this step if you want to use a single-domain, direct PCA basis. 
Currently we only support the skinning-correction basis for multi-domain
subspace simulation and subspace condensation simualtion. The results are
saved in the data path.

./bin/ComputeBasis ./cfg/chand.training.cfg
This call compute the subspace basis from the full space training data.
There are 3 parameters for computing the basis in the config file:

"partitioned basis": If 0, compute one basis for the entire mesh; If 1, compute
  a basis for each skeletal domain. The partitioning is done by
  assigning each vertex to the bone with maximum skinning weight.
"transformed basis": If 0, use the originals displacements for computing the
  basis; If 1, use the before-skinning displacement. As mentioned above, if you
  turn on "partitioned basis", "transformed basis" will be automatically set
  to 1 and an error message is emitted if TransformDisplacements hasn't been
  called yet.
"pca variance": The amount of variance you want to perserve when truncating
  the basis.

./bin/PartitionedInternalForceCubature ./cfg/chand.training.cfg
  This trains the material force cubatures for each skeletal domain. The single
  -domain version is InternalForceCubature.

./bin/SCFCubature ./cfg/chand.training.cfg
  This call trains the self-collision cubatures so that collisions already seen
  during the training can be quickly detected and resolved in subspace
  simulation.

./bin/PartitionedHybridSim ./cfg/chand.run.cfg
  This conducts subspace-condensation simulation on the hand calisthenic
  sequence, as in the SIGGRAPH 2015 video.

============================================
Documentation

The purpose of the current code release is to provide a working implementation
of the techniques from the above papers. While every effort has been made to 
make the code readable, it is still primarily a research prototype, and time 
constraints preclude its full documentation. Cubica++ follows the same coding
standard as Cubica (http://www.mat.ucsb.edu/~kim/cubica/docs.html) and the
code are well-organized into subdirectories based on their functionalities. 
Any of the main functions in the projects serves as a good entrance point to
understand the code. Timers (include/util/TIMING_BREAKDOWN.h) are inserted
throughout the code for profiling and most of the time the string inside
TIMING_BREAKDOWN::toc("...") explains what the chunk between tic() and toc()
is doing.
