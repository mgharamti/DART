#!/bin/bash
set -e

# Load the needed modules
module purge
module load ncarenv/23.09 nco/5.2.4 intel/2024.0.2 cray-mpich/8.1.27 \
            netcdf/4.9.2 conda/latest craype/2.7.23 ncview/2.1.9 \
            ncarcompilers/1.0.0 hdf5/1.14.3

# Activate the conda environment
conda activate /glade/work/gharamti/conda-envs/roms_dart_training

if ! jupyter kernelspec list | grep -q roms_dart_training ; then
    python -m ipykernel install --user \
        --name roms_dart_training \
        --display-name "ROMS-DART Training"
fi

# Test 
python - << EOF
import sys
import pydartdiags
print("Python:", sys.executable)
print("pydartdiags:", pydartdiags.__file__)
print("ROMS-DART Training kernel installed successfully")
EOF

jupyter kernelspec list

echo "ROMS-DART Training kernel is ready."
