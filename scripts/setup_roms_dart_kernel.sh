#!/bin/bash
set -e

# Activate the conda environment
conda activate /glade/work/gharamti/conda-envs/roms_dart_training

if ! jupyter kernelspec list | awk '{print $1}' | grep -qx "roms_dart_training" ; then
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

echo
jupyter kernelspec list

echo 
echo "ROMS-DART Training kernel is ready."
echo
