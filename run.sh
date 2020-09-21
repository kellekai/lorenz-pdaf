#!/bin/bash

#								 #
# INITIALIZATION #
#								 #

# generate observations
mpirun -n 4 ./a.out \
	-seed `date +"%N"` \
	-obs_share 5 \
	-obs_block 8 \
	-obs_gen 1 

# generate ensemble
for i in `seq 1 9`; do 
	mpirun -n 4 ./a.out \
		-seed `date +"%N"` \
		-member $i \
		-obs_share 5 \
		-obs_block 8 \
		-obs_gen 0 
done

#						 #
# assimilate #
#						 #

cd pdaf_from_scratch
mpirun -n 4 ./PDAF_offline
cd -

#								 		#
# GENERATE FORECAST #
#								 		#

# generate observations
mpirun -n 4 ./a.out \
	-seed `date +"%N"` \
	-obs_share 5 \
	-obs_block 8 \
	-obs_gen 1 \
  -epoch 1

# generate ensemble
for i in `seq 1 9`; do 
	mpirun -n 4 ./a.out \
		-seed `date +"%N"` \
		-member $i \
		-obs_share 5 \
		-obs_block 8 \
		-obs_gen 0 \
  	-epoch 1
done
