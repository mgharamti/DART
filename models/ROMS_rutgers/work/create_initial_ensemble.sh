#!/bin/bash

#PBS -N PSI
#PBS -e psi.err
#PBS -o psi.out
#PBS -l select=1:ncpus=128
#PBS -l walltime=00:10:00
#PBS -A NHAP0012 
#PBS -q main
#PBS -l job_priority=premium

work_dir=/glade/derecho/scratch/${USER}/inacawo/DART_training/models/ROMS_rutgers/work

# Cleanup
rm -f $work_dir/dart_log.out $work_dir/psi.* $work_dir/core 

# Run perturb_single_instance
mpirun -n 100 $work_dir/perturb_single_instance

# Move the ensemble 
mkdir -p ${work_dir}/small_ensemble
mv ${work_dir}/roms_mem*.nc ${work_dir}/small_ensemble
