#!/bin/bash
./bin/ConvertTetGen ./cfg/chand.fullsim.0.cfg
./bin/RigMesh ./cfg/chand.fullsim.0.cfg
./fullsim.sh
./bin/TransformDisplacements ./cfg/chand.training.cfg
./bin/ComputeBasis ./cfg/chand.training.cfg
./bin/PartitionedInternalForceCubature ./cfg/chand.training.cfg
./bin/SCFCubature ./cfg/chand.training.cfg
./bin/PartitionedHybridSim ./cfg/chand.run.cfg