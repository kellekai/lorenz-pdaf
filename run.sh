#!/bin/bash

APP=model
MAX_EPOCH=100
SHARE=20
BLOCK=4

#================#
# INITIALIZATION #
#================#

######################### 
# generate observations #
######################### 
mpirun -n 4 ./$APP \
	-seed `date +"%N"` \
	-obs_share $SHARE \
	-obs_block $BLOCK \
	-obs_gen T 

#####################
# generate ensemble #
#####################
for i in `seq 1 9`; do 
	mpirun -n 4 ./$APP \
		-seed `date +"%N"` \
		-member $i \
		-obs_share $SHARE \
		-obs_block $BLOCK \
		-obs_gen F 
done

##############
# ASSIMILATE #
##############
cd pdaf
mpirun -n 4 ./PDAF_offline \
  		-obs_share $SHARE \
  		-obs_block $BLOCK \
      -epoch 0
cd -

#====================#
# ASSIMILATION CYCLE #
#====================#

for epoch in `seq 1 $MAX_EPOCH`; do
  # generate observations
  mpirun -n 4 ./$APP \
  	-seed `date +"%N"` \
  	-obs_share $SHARE \
  	-obs_block $BLOCK \
  	-obs_gen T \
    -epoch $epoch
  
  # generate ensemble
  for i in `seq 1 9`; do 
  	mpirun -n 4 ./$APP \
  		-seed `date +"%N"` \
  		-member $i \
  		-obs_share $SHARE \
  		-obs_block $BLOCK \
  		-obs_gen F \
      -epoch $epoch
  done
  
  # assimilate
  cd pdaf
  mpirun -n 4 ./PDAF_offline \
  		-obs_share $SHARE \
  		-obs_block $BLOCK \
      -epoch $epoch
  cd -
done
