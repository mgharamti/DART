#!/bin/bash

# Load the needed modules
module purge
module load ncarenv/23.09 nco/5.2.4 intel/2024.0.2 cray-mpich/8.1.27 \
            netcdf/4.9.2 conda/latest craype/2.7.23 ncview/2.1.9 \
            ncarcompilers/1.0.0 hdf5/1.14.3

module list
