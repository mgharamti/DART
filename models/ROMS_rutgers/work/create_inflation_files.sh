#!/bin/bash

#PBS -N INF
#PBS -e inf.err
#PBS -o inf.out
#PBS -l select=1:ncpus=128
#PBS -l walltime=00:10:00
#PBS -A NHAP0012 
#PBS -q main
#PBS -l job_priority=premium

work_dir=/glade/derecho/scratch/${USER}/inacawo/DART_training/models/ROMS_rutgers/work

# Cleanup
rm -f $work_dir/dart_log.out $work_dir/inf.* $work_dir/core 

# Run fill_inflation_restart
cd $work_dir

./fill_inflation_restart
