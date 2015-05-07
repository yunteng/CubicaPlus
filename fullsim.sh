#!/bin/bash
for i in `seq 0 7`;
do
  ./bin/FullSim ./cfg/chand.fullsim.$i.cfg
done
