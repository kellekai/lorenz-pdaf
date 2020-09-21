#!/bin/bash

APP=model
MAX_EPOCH=2

#================#
# INITIALIZATION #
#================#

######################### 
# generate observations #
######################### 
mpirun -n 4 ./$APP \
	-seed `date +"%N"` \
	-obs_share 5 \
	-obs_block 8 \
	-obs_gen 1 

#####################
# generate ensemble #
#####################
for i in `seq 1 9`; do 
	mpirun -n 4 ./$APP \
		-seed `date +"%N"` \
		-member $i \
		-obs_share 5 \
		-obs_block 8 \
		-obs_gen 0 
done

##############
# ASSIMILATE #
##############
cd pdaf
mpirun -n 4 ./PDAF_offline \
  		-obs_share 5 \
  		-obs_block 8
cd -

#====================#
# ASSIMILATION CYCLE #
#====================#

for epoch in `seq 1 $MAX_EPOCH`; do
  # generate observations
  mpirun -n 4 ./$APP \
  	-seed `date +"%N"` \
  	-obs_share 5 \
  	-obs_block 8 \
  	-obs_gen 1 \
    -epoch $epoch
  
  # generate ensemble
  for i in `seq 1 9`; do 
  	mpirun -n 4 ./$APP \
  		-seed `date +"%N"` \
  		-member $i \
  		-obs_share 5 \
  		-obs_block 8 \
  		-obs_gen 0 \
      -epoch $epoch
  done
  
  # assimilate
  cd pdaf
  mpirun -n 4 ./PDAF_offline \
  		-obs_share 5 \
  		-obs_block 8
  cd -
done
