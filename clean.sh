#!bin/bash
export PDAF_ARCH=linux_gfortran_openmpi
cd pdaf
make cleanall
cd -
make cleanall
